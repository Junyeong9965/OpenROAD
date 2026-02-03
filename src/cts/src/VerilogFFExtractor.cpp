// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "VerilogFFExtractor.h"

#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>

#include "odb/db.h"
#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "sta/ExceptionPath.hh"
#include "sta/Graph.hh"
#include "sta/MinMax.hh"
#include "sta/Network.hh"
#include "sta/Path.hh"
#include "sta/PathEnd.hh"
#include "sta/PathExpanded.hh"
#include "sta/Sta.hh"
#include "utl/Logger.h"

namespace cts {

VerilogFFExtractor::VerilogFFExtractor(const std::string& verilog_path,
                                       odb::dbBlock* block,
                                       sta::dbSta* sta,
                                       sta::dbNetwork* network,
                                       utl::Logger* logger)
    : verilog_path_(verilog_path), block_(block), sta_(sta), network_(network), logger_(logger)
{
}

bool VerilogFFExtractor::isRegisterCell(const std::string& cell_type) const
{
  // Check if cell type starts with DFF or DFFH (common register prefixes)
  return cell_type.rfind("DFF", 0) == 0 || cell_type.rfind("DFFH", 0) == 0;
}

void VerilogFFExtractor::parseVerilog()
{
  logger_->info(utl::CTS, 350, "Parsing verilog file: {}", verilog_path_);

  std::ifstream file(verilog_path_);
  if (!file.is_open()) {
    logger_->error(utl::CTS, 351, "Cannot open verilog file: {}", verilog_path_);
  }

  // Read entire file
  std::stringstream buffer;
  buffer << file.rdbuf();
  std::string content = buffer.str();
  file.close();

  // Find module
  std::regex module_regex(R"(module\s+(\w+))");
  std::smatch module_match;
  if (!std::regex_search(content, module_match, module_regex)) {
    logger_->error(utl::CTS, 352, "No module found in verilog file");
  }

  std::string module_name = module_match[1].str();
  logger_->info(utl::CTS, 304, "Module: {}", module_name);

  // Find module body
  size_t start = module_match.position() + module_match.length();
  size_t endmod = content.find("endmodule", start);
  if (endmod == std::string::npos) {
    endmod = content.length();
  }
  std::string body = content.substr(start, endmod - start);

  // Parse all instances
  // Pattern: cell_type inst_name ( .pin(net), ... );
  // Handle escaped names like \dcnt[0]$_SDFFE_PN0P_ and hierarchical \u0.r0.out[24]$_SDFF_PP1_
  std::regex inst_regex(
      R"((\w+)\s+(\\?[\w\.\[\]$_]+)\s*\(([^;]+)\);)",
      std::regex_constants::multiline);

  auto inst_begin = std::sregex_iterator(body.begin(), body.end(), inst_regex);
  auto inst_end = std::sregex_iterator();

  for (std::sregex_iterator i = inst_begin; i != inst_end; ++i) {
    std::smatch match = *i;
    std::string cell_type = match[1].str();
    std::string inst_name = match[2].str();
    std::string conn_block = match[3].str();

    // Trim whitespace from inst_name
    inst_name.erase(0, inst_name.find_first_not_of(" \t\n\r"));
    inst_name.erase(inst_name.find_last_not_of(" \t\n\r") + 1);

    // Parse pin connections: .pin_name(net_name)
    std::regex pin_regex(R"(\.(\w+)\s*\(\s*([^\)]+)\s*\))");
    auto pin_begin = std::sregex_iterator(conn_block.begin(), conn_block.end(), pin_regex);
    auto pin_end = std::sregex_iterator();

    RegisterInfo reg_info;
    reg_info.cell_type = cell_type;

    for (std::sregex_iterator j = pin_begin; j != pin_end; ++j) {
      std::smatch pin_match = *j;
      std::string pin_name = pin_match[1].str();
      std::string net_name = pin_match[2].str();

      // Trim whitespace from net_name
      net_name.erase(0, net_name.find_first_not_of(" \t\n\r"));
      net_name.erase(net_name.find_last_not_of(" \t\n\r") + 1);

      reg_info.pins[pin_name] = net_name;
      net_connections_[net_name].push_back({inst_name, pin_name});
    }

    // Store all instances
    all_instances_[inst_name] = reg_info;

    // Also store in registers_ if it's a register
    if (isRegisterCell(cell_type)) {
      registers_[inst_name] = reg_info;
    }
  }

  logger_->info(utl::CTS, 305, "Total registers found: {}", registers_.size());
  logger_->info(utl::CTS, 306, "Total nets: {}", net_connections_.size());
  logger_->info(utl::CTS, 321, "Total instances parsed: {}", all_instances_.size());
}

std::vector<FFEdgeVerilog> VerilogFFExtractor::extractFFEdges()
{
  logger_->info(utl::CTS, 307, "Extracting FF-to-FF edges...");

  std::vector<FFEdgeVerilog> edges;

  // Debug: print first few registers
  int debug_count = 0;
  for (const auto& [inst_name, info] : registers_) {
    if (debug_count < 5) {
      std::string pins_str;
      for (const auto& [pin, net] : info.pins) {
        pins_str += pin + "->" + net + " ";
      }
      logger_->info(utl::CTS, 320, "Register {}: {}", inst_name, pins_str);
      debug_count++;
    }
  }

  // For each register, find its Q output net
  for (const auto& [from_inst, from_info] : registers_) {
    const auto& pins = from_info.pins;

    // Find Q or QN output
    std::string q_net;
    auto q_it = pins.find("Q");
    if (q_it != pins.end()) {
      q_net = q_it->second;
    } else {
      auto qn_it = pins.find("QN");
      if (qn_it != pins.end()) {
        q_net = qn_it->second;
      }
    }

    if (q_net.empty()) {
      continue;  // Some FFs may have only complemented outputs
    }

    // BFS to find all reachable FF D pins through combinational logic
    std::unordered_set<std::string> visited_nets;
    std::vector<std::string> queue;
    queue.push_back(q_net);
    visited_nets.insert(q_net);

    while (!queue.empty()) {
      std::string current_net = queue.back();
      queue.pop_back();

      auto net_it = net_connections_.find(current_net);
      if (net_it == net_connections_.end()) {
        continue;
      }

      const auto& connected = net_it->second;

      for (const auto& [to_inst, to_pin] : connected) {
        if (to_inst == from_inst) {
          continue;  // Skip self-loops
        }

        // Check if destination is a register D pin
        if (registers_.find(to_inst) != registers_.end()) {
          if (to_pin == "D") {
            FFEdgeVerilog edge;
            edge.from_ff = from_inst;
            edge.to_ff = to_inst;
            edge.via_net = current_net;

            // Fill locations from ODB if available
            fillLocations(edge);

            edges.push_back(edge);
          }
          // Don't traverse through registers
          continue;
        }

        // For non-register instances (combinational gates), find their outputs
        auto gate_it = all_instances_.find(to_inst);
        if (gate_it != all_instances_.end()) {
          // Find output pins: Y, Z, Q, QN, ZN (standard), CON, SN (adder outputs)
          // Note: 'C' is NOT an output - it's an input pin in OA211, AO211, etc.
          const auto& gate_pins = gate_it->second.pins;
          for (const auto& [gate_pin, gate_net] : gate_pins) {
            // Common output pin names (standard gates + adders)
            if (gate_pin == "Y" || gate_pin == "Z" || gate_pin == "Q" ||
                gate_pin == "QN" || gate_pin == "ZN" ||
                gate_pin == "CON" || gate_pin == "SN") {
              if (visited_nets.find(gate_net) == visited_nets.end()) {
                queue.push_back(gate_net);
                visited_nets.insert(gate_net);
              }
            }
          }
        }
      }
    }
  }

  logger_->info(utl::CTS, 308, "Found {} FF-to-FF edges", edges.size());

  return edges;
}

// Helper function to escape brackets for ODB lookup
static std::string escapeForODB(const std::string& name)
{
  std::string result;
  std::string input = name;

  // Remove leading backslash if present
  if (!input.empty() && input[0] == '\\') {
    input = input.substr(1);
  }

  // Escape brackets: [ -> \[, ] -> \]
  for (char c : input) {
    if (c == '[' || c == ']') {
      result += '\\';
    }
    result += c;
  }
  return result;
}

void VerilogFFExtractor::fillLocations(FFEdgeVerilog& edge)
{
  if (!block_) {
    return;
  }

  // Get locations from ODB (need to escape brackets)
  std::string from_odb = escapeForODB(edge.from_ff);
  std::string to_odb = escapeForODB(edge.to_ff);

  odb::dbInst* from_inst = block_->findInst(from_odb.c_str());
  odb::dbInst* to_inst = block_->findInst(to_odb.c_str());

  if (from_inst) {
    int x, y;
    from_inst->getLocation(x, y);
    edge.from_x = x;
    edge.from_y = y;
    // Get tier from master name (_upper -> 1, _bottom -> 0)
    std::string master = from_inst->getMaster()->getName();
    edge.from_tier = (master.find("_upper") != std::string::npos) ? 1 : 0;
  }

  if (to_inst) {
    int x, y;
    to_inst->getLocation(x, y);
    edge.to_x = x;
    edge.to_y = y;
    // Get tier from master name (_upper -> 1, _bottom -> 0)
    std::string master = to_inst->getMaster()->getName();
    edge.to_tier = (master.find("_upper") != std::string::npos) ? 1 : 0;
  }
}

void VerilogFFExtractor::fillTimingInfo(std::vector<FFEdgeVerilog>& edges)
{
  if (!sta_ || !network_) {
    logger_->warn(utl::CTS, 322, "STA not available, skipping timing info");
    return;
  }

  logger_->info(utl::CTS, 323, "Filling timing info for {} edges...", edges.size());

  // Ensure timing is updated
  sta_->ensureGraph();
  sta_->ensureClkArrivals();
  sta_->updateTiming(false);

  sta::Graph* graph = sta_->graph();
  if (!graph) {
    logger_->warn(utl::CTS, 324, "STA graph not available");
    return;
  }

  constexpr float SEC_TO_NS = 1e9f;
  int timing_found = 0;
  int hold_found = 0;

  // Group edges by from_ff for efficient findPathEnds calls
  std::unordered_map<std::string, std::vector<size_t>> from_ff_edges;
  for (size_t i = 0; i < edges.size(); i++) {
    from_ff_edges[edges[i].from_ff].push_back(i);
  }

  logger_->info(utl::CTS, 327, "Grouped into {} unique from_ff startpoints", from_ff_edges.size());

  // Build edge lookup: (from_ff, to_ff) -> edge index
  std::unordered_map<std::string, size_t> edge_lookup;
  for (size_t i = 0; i < edges.size(); i++) {
    std::string key = edges[i].from_ff + "|" + edges[i].to_ff;
    edge_lookup[key] = i;
  }

  int call_count = 0;

  // Helper lambda to find edge index from ODB instance name
  auto findEdgeIdx = [&](const std::string& from_ff_name,
                         const std::string& odb_name,
                         const std::vector<size_t>& edge_indices) -> int {
    std::string to_name;
    std::string key;
    auto it = edge_lookup.end();

    // Format 1: \name (Verilog escaped)
    to_name = "\\" + odb_name;
    key = from_ff_name + "|" + to_name;
    it = edge_lookup.find(key);

    // Format 2: name (raw ODB)
    if (it == edge_lookup.end()) {
      to_name = odb_name;
      key = from_ff_name + "|" + to_name;
      it = edge_lookup.find(key);
    }

    // Format 3: Try matching with escapeForODB conversion
    if (it == edge_lookup.end()) {
      for (size_t idx : edge_indices) {
        std::string edge_to_odb = escapeForODB(edges[idx].to_ff);
        if (edge_to_odb == odb_name) {
          it = edge_lookup.find(from_ff_name + "|" + edges[idx].to_ff);
          break;
        }
      }
    }

    return (it != edge_lookup.end()) ? static_cast<int>(it->second) : -1;
  };

  // For each unique from_ff, call findPathEnds for both SETUP and HOLD
  for (const auto& [from_ff_name, edge_indices] : from_ff_edges) {
    // Find from_ff instance in ODB
    std::string from_odb_name = escapeForODB(from_ff_name);
    odb::dbInst* from_odb_inst = block_ ? block_->findInst(from_odb_name.c_str()) : nullptr;
    if (!from_odb_inst) continue;

    sta::Instance* from_sta_inst = network_->dbToSta(from_odb_inst);
    if (!from_sta_inst) continue;

    // ============ SETUP (max) timing ============
    sta::InstanceSet* from_insts_setup = new sta::InstanceSet(network_);
    from_insts_setup->insert(from_sta_inst);
    sta::ExceptionFrom* from_setup = sta_->makeExceptionFrom(
        nullptr, nullptr, from_insts_setup, sta::RiseFallBoth::riseFall());

    sta::PathEndSeq path_ends_setup = sta_->findPathEnds(
        from_setup, nullptr, nullptr, true, nullptr,
        sta::MinMaxAll::max(), 10000, 1, true, false,
        -sta::INF, sta::INF, false, nullptr,
        true, false, false, false, false, false);

    // Process SETUP path ends
    for (sta::PathEnd* path_end : path_ends_setup) {
      if (!path_end) continue;
      const sta::Path* data_path = path_end->path();
      if (!data_path) continue;

      sta::Pin* end_pin = data_path->pin(sta_);
      if (!end_pin) continue;

      sta::Instance* end_inst = network_->instance(end_pin);
      if (!end_inst) continue;

      odb::dbInst* end_odb_inst = network_->staToDb(end_inst);
      if (!end_odb_inst) continue;

      std::string odb_name = end_odb_inst->getName();
      int idx = findEdgeIdx(from_ff_name, odb_name, edge_indices);

      if (idx >= 0 && !edges[idx].has_timing) {
        edges[idx].slack_max = sta::delayAsFloat(path_end->slack(sta_)) * SEC_TO_NS;
        edges[idx].arrival_max = sta::delayAsFloat(data_path->arrival()) * SEC_TO_NS;
        edges[idx].required_max = sta::delayAsFloat(path_end->requiredTime(sta_)) * SEC_TO_NS;
        edges[idx].has_timing = true;
        timing_found++;
      }
    }

    // ============ HOLD (min) timing ============
    sta::InstanceSet* from_insts_hold = new sta::InstanceSet(network_);
    from_insts_hold->insert(from_sta_inst);
    sta::ExceptionFrom* from_hold = sta_->makeExceptionFrom(
        nullptr, nullptr, from_insts_hold, sta::RiseFallBoth::riseFall());

    sta::PathEndSeq path_ends_hold = sta_->findPathEnds(
        from_hold, nullptr, nullptr, true, nullptr,
        sta::MinMaxAll::min(), 10000, 1, true, false,
        -sta::INF, sta::INF, false, nullptr,
        false, true, false, false, false, false);

    // Process HOLD path ends - use vertex->pin() instead of path->pin()
    for (sta::PathEnd* path_end : path_ends_hold) {
      if (!path_end) continue;

      // Get endpoint vertex directly from PathEnd
      sta::Vertex* end_vertex = path_end->vertex(sta_);
      if (!end_vertex) continue;

      sta::Pin* end_pin = end_vertex->pin();
      if (!end_pin) continue;

      sta::Instance* end_inst = network_->instance(end_pin);
      if (!end_inst) continue;

      odb::dbInst* end_odb_inst = network_->staToDb(end_inst);
      if (!end_odb_inst) continue;

      std::string odb_name = end_odb_inst->getName();
      int idx = findEdgeIdx(from_ff_name, odb_name, edge_indices);

      if (idx >= 0) {
        const sta::Path* data_path = path_end->path();
        if (data_path) {
          edges[idx].slack_min = sta::delayAsFloat(path_end->slack(sta_)) * SEC_TO_NS;
          edges[idx].arrival_min = sta::delayAsFloat(data_path->arrival()) * SEC_TO_NS;
          edges[idx].required_min = sta::delayAsFloat(path_end->requiredTime(sta_)) * SEC_TO_NS;
          hold_found++;
        }
      }
    }

    call_count++;
  }

  logger_->info(utl::CTS, 325, "Found setup timing for {} / {} edges, hold timing for {} edges ({} calls)",
                timing_found, edges.size(), hold_found, call_count);
}

void VerilogFFExtractor::writeCSV(const std::vector<FFEdgeVerilog>& edges,
                                  const std::string& output_file)
{
  std::ofstream file(output_file);
  if (!file.is_open()) {
    logger_->error(utl::CTS, 309, "Cannot open output file: {}", output_file);
  }

  // Write header (with setup and hold timing + tier info)
  file << "from_ff,to_ff,slack_max_ns,slack_min_ns,arrival_max_ns,arrival_min_ns,required_max_ns,required_min_ns,from_x,from_y,from_tier,to_x,to_y,to_tier\n";

  // Write edges
  for (const auto& edge : edges) {
    file << edge.from_ff << "," << edge.to_ff << ","
         << edge.slack_max << "," << edge.slack_min << ","
         << edge.arrival_max << "," << edge.arrival_min << ","
         << edge.required_max << "," << edge.required_min << ","
         << edge.from_x << "," << edge.from_y << "," << edge.from_tier << ","
         << edge.to_x << "," << edge.to_y << "," << edge.to_tier << "\n";
  }

  file.close();

  logger_->info(
      utl::CTS, 310, "Wrote {} edges to {}", edges.size(), output_file);
}

}  // namespace cts
