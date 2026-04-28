// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "VerilogFFExtractor.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <regex>
#include <sstream>

#include "Cts3DDatabase.h"  // Added for tier queries
#include "odb/db.h"
#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "sta/ExceptionPath.hh"
#include "sta/Graph.hh"
#include "sta/Liberty.hh"   // BUF_MACRO: for hasSequentials()
#include "sta/MinMax.hh"
#include "sta/Network.hh"
#include "sta/Path.hh"
#include "sta/PathEnd.hh"
#include "sta/PathExpanded.hh"
#include "sta/Sta.hh"
#include "utl/Logger.h"

namespace cts {

namespace {

bool isNamedHardMacroMaster(const odb::dbMaster* master)
{
  if (master == nullptr) {
    return false;
  }
  std::string name = master->getName();
  std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) {
    return std::tolower(c);
  });
  return name.rfind("fakeram", 0) == 0 || name.rfind("fakeregfile", 0) == 0
         || name.rfind("sram_", 0) == 0;
}

bool isMacroBlockInst(const odb::dbInst* inst)
{
  if (inst == nullptr) {
    return false;
  }
  if (inst->isBlock()) {
    return true;
  }
  const odb::dbMaster* master = inst->getMaster();
  return (master != nullptr && master->isBlock())
         || isNamedHardMacroMaster(master);
}

}  // namespace

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
  logger_->info(utl::CTS, 354, "Parsing verilog file: {}", verilog_path_);

  std::ifstream file(verilog_path_);
  if (!file.is_open()) {
    logger_->error(utl::CTS, 355, "Cannot open verilog file: {}", verilog_path_);
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
    logger_->error(utl::CTS, 356, "No module found in verilog file");
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
    // Replaced master name parsing with Cts3DDatabase query
    // Fix fallback: master names use "_upper" not "__upper"
    if (cts3dDb_) {
      edge.from_tier = cts3dDb_->getInstTier(from_inst);
    } else {
      std::string master = from_inst->getMaster()->getName();
      edge.from_tier = (master.find("_upper") != std::string::npos) ? 1 : 0;
    }
  }

  if (to_inst) {
    int x, y;
    to_inst->getLocation(x, y);
    edge.to_x = x;
    edge.to_y = y;
    // Replaced master name parsing with Cts3DDatabase query
    // Fix fallback: master names use "_upper" not "__upper"
    if (cts3dDb_) {
      edge.to_tier = cts3dDb_->getInstTier(to_inst);
    } else {
      std::string master = to_inst->getMaster()->getName();
      edge.to_tier = (master.find("_upper") != std::string::npos) ? 1 : 0;
    }
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

  // Read PATH_LIMIT env var (same as fillTimingInfoBatch)
  int n_path_limit = 10000;
  const char* plimit_env_ti = std::getenv("PATH_LIMIT");
  if (plimit_env_ti) {
    n_path_limit = std::max(10000, std::atoi(plimit_env_ti));
  }

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
        sta::MinMaxAll::max(), n_path_limit, 1, true, false,
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
        sta::MinMaxAll::min(), n_path_limit, 1, true, false,
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

// Source-side batched timing fill.
// Same logic as fillTimingInfo but batches K source FFs into one InstanceSet
// per findPathEnds call. Reduces per-call overhead (graph setup, exception
// creation) while preserving per-edge exact results.
// BATCH_SIZE env var controls K (default 10).
// PARTITION_FILE env var optionally filters source FFs for parallel mode.
void VerilogFFExtractor::fillTimingInfoBatch(std::vector<FFEdgeVerilog>& edges)
{
  if (!sta_ || !network_) {
    logger_->warn(utl::CTS, 640, "STA not available, skipping batch timing fill");
    return;
  }

  int batch_size = 10;
  const char* batch_env = std::getenv("BATCH_SIZE");
  if (batch_env) {
    batch_size = std::max(1, std::atoi(batch_env));
  }

  logger_->info(utl::CTS, 641,
      "Batch timing fill for {} edges (source-side, batch_size={})",
      edges.size(), batch_size);

  sta_->ensureGraph();
  sta_->ensureClkArrivals();
  sta_->updateTiming(false);

  constexpr float SEC_TO_NS = 1e9f;
  int timing_found = 0;
  int hold_found = 0;
  int call_count = 0;

  // Group edges by source FF (from_ff) — same as fillTimingInfo
  std::unordered_map<std::string, std::vector<size_t>> from_ff_edges;
  for (size_t i = 0; i < edges.size(); i++) {
    from_ff_edges[edges[i].from_ff].push_back(i);
  }

  // Optional partition filter
  std::unordered_set<std::string> partition_ffs;
  bool use_partition = false;
  const char* part_file_env = std::getenv("PARTITION_FILE");
  if (part_file_env && std::string(part_file_env) != "") {
    std::ifstream part_file(part_file_env);
    if (part_file.is_open()) {
      std::string line;
      while (std::getline(part_file, line)) {
        if (!line.empty()) partition_ffs.insert(line);
      }
      use_partition = true;
      logger_->info(utl::CTS, 648, "Partition filter: {} source FFs from {}",
                    partition_ffs.size(), part_file_env);
    }
  }

  // Build edge lookup: (from_ff, to_ff) -> edge index
  std::unordered_map<std::string, size_t> edge_lookup;
  for (size_t i = 0; i < edges.size(); i++) {
    std::string key = edges[i].from_ff + "|" + edges[i].to_ff;
    edge_lookup[key] = i;
  }

  // Resolve source FFs to STA instances, optionally filtered by partition
  struct SourceFF {
    std::string name;
    sta::Instance* sta_inst;
  };
  std::vector<SourceFF> source_ffs;
  for (const auto& [ff_name, _] : from_ff_edges) {
    if (use_partition && partition_ffs.find(ff_name) == partition_ffs.end()) {
      continue;
    }
    std::string odb_name = escapeForODB(ff_name);
    odb::dbInst* odb_inst = block_ ? block_->findInst(odb_name.c_str()) : nullptr;
    if (!odb_inst) continue;
    sta::Instance* sta_inst = network_->dbToSta(odb_inst);
    if (!sta_inst) continue;
    source_ffs.push_back({ff_name, sta_inst});
  }

  logger_->info(utl::CTS, 643,
      "Resolved {} source FFs ({} total unique, {} after partition filter)",
      source_ffs.size(), from_ff_edges.size(),
      use_partition ? partition_ffs.size() : from_ff_edges.size());

  // Reuse findEdgeIdx from fillTimingInfo — match endpoint ODB name to edge
  auto findEdgeIdx = [&](const std::string& from_ff_name,
                         const std::string& odb_name,
                         const std::vector<size_t>& edge_indices) -> int {
    std::string key;

    // Format 1: escaped Verilog name
    key = from_ff_name + "|" + "\\" + odb_name;
    auto it = edge_lookup.find(key);
    if (it != edge_lookup.end()) return static_cast<int>(it->second);

    // Format 2: raw ODB name
    key = from_ff_name + "|" + odb_name;
    it = edge_lookup.find(key);
    if (it != edge_lookup.end()) return static_cast<int>(it->second);

    // Format 3: escape ODB name to match edge to_ff
    for (size_t idx : edge_indices) {
      std::string edge_to_odb = escapeForODB(edges[idx].to_ff);
      if (edge_to_odb == odb_name) {
        return static_cast<int>(idx);
      }
    }
    return -1;
  };

  int n_path_limit = 10000;
  const char* plimit_env = std::getenv("PATH_LIMIT");
  if (plimit_env) {
    n_path_limit = std::max(10000, std::atoi(plimit_env));
  }

  // Process source FFs in batches of batch_size
  for (size_t batch_start = 0; batch_start < source_ffs.size();
       batch_start += batch_size)
  {
    size_t batch_end = std::min(batch_start + batch_size, source_ffs.size());

    // ============ SETUP (max) timing ============
    {
      sta::InstanceSet* from_insts = new sta::InstanceSet(network_);
      for (size_t b = batch_start; b < batch_end; b++) {
        from_insts->insert(source_ffs[b].sta_inst);
      }
      sta::ExceptionFrom* from = sta_->makeExceptionFrom(
          nullptr, nullptr, from_insts, sta::RiseFallBoth::riseFall());

      sta::PathEndSeq path_ends = sta_->findPathEnds(
          from, nullptr, nullptr, true, nullptr,
          sta::MinMaxAll::max(), n_path_limit, 1, true, false,
          -sta::INF, sta::INF, false, nullptr,
          true, false, false, false, false, false);

      if (static_cast<int>(path_ends.size()) >= n_path_limit) {
        logger_->warn(utl::CTS, 649,
            "Setup findPathEnds hit path limit {} — increase PATH_LIMIT env var",
            n_path_limit);
      }

      // For each returned path, find which source FF it belongs to
      for (sta::PathEnd* pe : path_ends) {
        if (!pe) continue;
        const sta::Path* data_path = pe->path();
        if (!data_path) continue;

        // Get capture FF (endpoint)
        sta::Pin* end_pin = data_path->pin(sta_);
        if (!end_pin) continue;
        sta::Instance* end_inst = network_->instance(end_pin);
        if (!end_inst) continue;
        odb::dbInst* end_odb = network_->staToDb(end_inst);
        if (!end_odb) continue;
        std::string capture_odb = end_odb->getName();

        // Get launch FF from path expansion
        sta::PathExpanded expanded(data_path, sta_);
        if (expanded.size() == 0) continue;
        const sta::Path* start_path = expanded.startPath();
        if (!start_path) start_path = expanded.path(0);
        if (!start_path) continue;
        sta::Pin* start_pin = start_path->pin(sta_);
        if (!start_pin || network_->isTopLevelPort(start_pin)) continue;
        sta::Instance* start_inst = network_->instance(start_pin);
        if (!start_inst) continue;
        odb::dbInst* start_odb = network_->staToDb(start_inst);
        if (!start_odb) continue;
        std::string launch_odb = start_odb->getName();

        // Match to edge: try each source FF in batch
        for (size_t b = batch_start; b < batch_end; b++) {
          const auto& src = source_ffs[b];
          // Check if this path's launch FF matches this batch source FF
          std::string src_odb = escapeForODB(src.name);
          odb::dbInst* src_inst = block_->findInst(src_odb.c_str());
          if (!src_inst || std::string(src_inst->getName()) != launch_odb) continue;

          auto edge_it = from_ff_edges.find(src.name);
          if (edge_it == from_ff_edges.end()) continue;

          int idx = findEdgeIdx(src.name, capture_odb, edge_it->second);
          if (idx >= 0 && !edges[idx].has_timing) {
            edges[idx].slack_max = sta::delayAsFloat(pe->slack(sta_)) * SEC_TO_NS;
            edges[idx].arrival_max = sta::delayAsFloat(data_path->arrival()) * SEC_TO_NS;
            edges[idx].required_max = sta::delayAsFloat(pe->requiredTime(sta_)) * SEC_TO_NS;
            edges[idx].has_timing = true;
            timing_found++;
          }
          break;
        }
      }
      call_count++;
    }

    // ============ HOLD (min) timing ============
    {
      sta::InstanceSet* from_insts = new sta::InstanceSet(network_);
      for (size_t b = batch_start; b < batch_end; b++) {
        from_insts->insert(source_ffs[b].sta_inst);
      }
      sta::ExceptionFrom* from = sta_->makeExceptionFrom(
          nullptr, nullptr, from_insts, sta::RiseFallBoth::riseFall());

      sta::PathEndSeq path_ends = sta_->findPathEnds(
          from, nullptr, nullptr, true, nullptr,
          sta::MinMaxAll::min(), n_path_limit, 1, true, false,
          -sta::INF, sta::INF, false, nullptr,
          false, true, false, false, false, false);

      for (sta::PathEnd* pe : path_ends) {
        if (!pe) continue;

        sta::Vertex* end_vertex = pe->vertex(sta_);
        if (!end_vertex) continue;
        sta::Pin* end_pin = end_vertex->pin();
        if (!end_pin) continue;
        sta::Instance* end_inst = network_->instance(end_pin);
        if (!end_inst) continue;
        odb::dbInst* end_odb = network_->staToDb(end_inst);
        if (!end_odb) continue;
        std::string capture_odb = end_odb->getName();

        // Get launch FF
        const sta::Path* data_path = pe->path();
        if (!data_path) continue;
        sta::PathExpanded expanded(data_path, sta_);
        if (expanded.size() == 0) continue;
        const sta::Path* start_path = expanded.startPath();
        if (!start_path) start_path = expanded.path(0);
        if (!start_path) continue;
        sta::Pin* start_pin = start_path->pin(sta_);
        if (!start_pin || network_->isTopLevelPort(start_pin)) continue;
        sta::Instance* start_inst = network_->instance(start_pin);
        if (!start_inst) continue;
        odb::dbInst* start_odb = network_->staToDb(start_inst);
        if (!start_odb) continue;
        std::string launch_odb = start_odb->getName();

        for (size_t b = batch_start; b < batch_end; b++) {
          const auto& src = source_ffs[b];
          std::string src_odb = escapeForODB(src.name);
          odb::dbInst* src_inst = block_->findInst(src_odb.c_str());
          if (!src_inst || std::string(src_inst->getName()) != launch_odb) continue;

          auto edge_it = from_ff_edges.find(src.name);
          if (edge_it == from_ff_edges.end()) continue;

          int idx = findEdgeIdx(src.name, capture_odb, edge_it->second);
          if (idx >= 0 && !edges[idx].has_hold) {
            edges[idx].slack_min = sta::delayAsFloat(pe->slack(sta_)) * SEC_TO_NS;
            edges[idx].arrival_min = sta::delayAsFloat(data_path->arrival()) * SEC_TO_NS;
            edges[idx].required_min = sta::delayAsFloat(pe->requiredTime(sta_)) * SEC_TO_NS;
            edges[idx].has_hold = true;
            hold_found++;
          }
          break;
        }
      }
      call_count++;
    }

    // Progress report every 100 batches
    if ((batch_start / batch_size) % 100 == 0) {
      logger_->info(utl::CTS, 645,
          "  Batch progress: {} / {} source FFs, {} setup + {} hold found "
          "({} findPathEnds calls)",
          batch_end, source_ffs.size(), timing_found, hold_found, call_count);
    }
  }

  logger_->info(utl::CTS, 646,
      "Batch timing complete: {} setup, {} hold for {} edges "
      "({} findPathEnds calls, batch_size={})",
      timing_found, hold_found, edges.size(), call_count, batch_size);
}

void VerilogFFExtractor::writeCSV(const std::vector<FFEdgeVerilog>& edges,
                                  const std::string& output_file)
{
  std::ofstream file(output_file);
  if (!file.is_open()) {
    logger_->error(utl::CTS, 309, "Cannot open output file: {}", output_file);
  }

  // Write header (with setup and hold timing + tier info + clock net)
  file << "from_ff,to_ff,slack_max_ns,slack_min_ns,arrival_max_ns,arrival_min_ns,required_max_ns,required_min_ns,from_x,from_y,from_tier,to_x,to_y,to_tier,from_clock_net,to_clock_net\n";

  // Write edges
  for (const auto& edge : edges) {
    file << edge.from_ff << "," << edge.to_ff << ","
         << edge.slack_max << "," << edge.slack_min << ","
         << edge.arrival_max << "," << edge.arrival_min << ","
         << edge.required_max << "," << edge.required_min << ","
         << edge.from_x << "," << edge.from_y << "," << edge.from_tier << ","
         << edge.to_x << "," << edge.to_y << "," << edge.to_tier << ","
         << edge.from_clock_net << "," << edge.to_clock_net << "\n";
  }

  file.close();

  logger_->info(
      utl::CTS, 310, "Wrote {} edges to {}", edges.size(), output_file);
}

// ============================================================================
// Exhaustive IO timing edge extraction via STA.
//
// Replaces the Tcl extract_io_timing_edges_phase1a (worst-N truncation) with
// per-port findPathEnds calls.  Every (port, FF) pair with a valid STA timing
// path is captured — no group_path_count limit.
//
// PI→FF: for each input BTerm, findPathEnds (setup + hold).
// FF→PO: for each output BTerm, findPathEnds (setup + hold).
// ============================================================================
std::vector<IOEdgeVerilog> VerilogFFExtractor::extractIOEdges()
{
  std::vector<IOEdgeVerilog> result;

  if (!sta_ || !network_ || !block_) {
    logger_->warn(utl::CTS, 600,
                  "3D-CTS IO: STA/network/block not available, skipping IO extraction");
    return result;
  }

  sta_->ensureGraph();
  sta_->ensureClkArrivals();
  sta_->updateTiming(false);

  constexpr float SEC_TO_NS = 1e9f;

  // Build reverse lookup: ODB inst* → name (for FF name in IO CSV).
  // BUF_MACRO: use odb_registers_ if available (ODB-based),
  // fall back to Verilog-parsed registers_ for backward compatibility.
  std::unordered_map<odb::dbInst*, std::string> inst_to_name;
  if (!odb_registers_.empty()) {
    // ODB-based: odb_registers_ already has dbInst* → name mapping
    for (const auto& [name, inst] : odb_registers_) {
      inst_to_name[inst] = name;
    }
  } else {
    // Verilog-based fallback (legacy path)
    for (const auto& [verilog_name, info] : registers_) {
      std::string odb_name = escapeForODB(verilog_name);
      odb::dbInst* inst = block_->findInst(odb_name.c_str());
      if (inst) {
        inst_to_name[inst] = inst->getName();
      }
    }
  }

  // Collect input / output BTerms
  std::vector<odb::dbBTerm*> input_ports, output_ports;
  for (auto* bterm : block_->getBTerms()) {
    odb::dbSigType sig = bterm->getSigType();
    if (sig == odb::dbSigType::CLOCK) {
      continue;  // Skip clock ports — not data IO
    }
    auto io = bterm->getIoType();
    if (io == odb::dbIoType::INPUT) {
      input_ports.push_back(bterm);
    } else if (io == odb::dbIoType::OUTPUT) {
      output_ports.push_back(bterm);
    }
  }

  logger_->info(utl::CTS, 601,
                "3D-CTS IO: extracting IO edges ({} input ports, {} output ports)",
                input_ports.size(), output_ports.size());

  // Edge map: key = "type|port|ff" → IOEdgeVerilog (keep worst slack per pair)
  std::unordered_map<std::string, IOEdgeVerilog> edge_map;

  // Helper: get FF name from STA endpoint instance
  auto getFFName = [&](sta::Pin* pin) -> std::string {
    if (!pin) return "";
    sta::Instance* inst = network_->instance(pin);
    if (!inst) return "";
    odb::dbInst* odb_inst = network_->staToDb(inst);
    if (!odb_inst) return "";
    // Use inst_to_name if available, else ODB getName()
    auto it = inst_to_name.find(odb_inst);
    if (it != inst_to_name.end()) return it->second;
    return odb_inst->getName();
  };

  // Per-clock-net LP: build FF -> clock net cache for IO edges.
  std::unordered_map<std::string, std::string> io_ff_clock_cache;
  for (const auto& [name, inst] : odb_registers_) {
    io_ff_clock_cache[name] = traceClockNetRoot(inst);
  }

  int pi_setup_ct = 0, pi_hold_ct = 0, po_setup_ct = 0, po_hold_ct = 0;
  int pi_removal_ct = 0, pi_recovery_ct = 0;

  // ==================== PI→FF ====================
  for (auto* bterm : input_ports) {
    std::string port_name = bterm->getName();
    sta::Pin* port_pin = network_->dbToSta(bterm);
    if (!port_pin) continue;

    // --- SETUP (max) ---
    {
      sta::PinSet* from_pins = new sta::PinSet(network_);
      from_pins->insert(port_pin);
      sta::ExceptionFrom* from = sta_->makeExceptionFrom(
          from_pins, nullptr, nullptr, sta::RiseFallBoth::riseFall());

      sta::PathEndSeq path_ends = sta_->findPathEnds(
          from, nullptr, nullptr, true, nullptr,
          sta::MinMaxAll::max(), 10000, 1, true, false,
          -sta::INF, sta::INF, false, nullptr,
          true, false, false, false, false, false);

      for (sta::PathEnd* pe : path_ends) {
        if (!pe) continue;
        const sta::Path* path = pe->path();
        if (!path) continue;
        sta::Pin* end_pin = path->pin(sta_);
        std::string ff = getFFName(end_pin);
        if (ff.empty()) continue;

        std::string key = "PI_TO_FF|" + port_name + "|" + ff;
        auto& e = edge_map[key];
        e.edge_type = "PI_TO_FF";
        e.port_name = port_name;
        e.ff_name = ff;
        e.ff_clock_net = io_ff_clock_cache.count(ff) ? io_ff_clock_cache[ff] : "";
        float slack = sta::delayAsFloat(pe->slack(sta_)) * SEC_TO_NS;
        if (!e.has_setup || slack < e.slack_setup_ns) {
          e.slack_setup_ns = slack;
          e.has_setup = true;
          pi_setup_ct++;
        }
      }
    }

    // --- HOLD (min) ---
    {
      sta::PinSet* from_pins = new sta::PinSet(network_);
      from_pins->insert(port_pin);
      sta::ExceptionFrom* from = sta_->makeExceptionFrom(
          from_pins, nullptr, nullptr, sta::RiseFallBoth::riseFall());

      sta::PathEndSeq path_ends = sta_->findPathEnds(
          from, nullptr, nullptr, true, nullptr,
          sta::MinMaxAll::min(), 10000, 1, true, false,
          -sta::INF, sta::INF, false, nullptr,
          false, true, false, false, false, false);

      for (sta::PathEnd* pe : path_ends) {
        if (!pe) continue;
        sta::Vertex* end_vertex = pe->vertex(sta_);
        if (!end_vertex) continue;
        sta::Pin* end_pin = end_vertex->pin();
        std::string ff = getFFName(end_pin);
        if (ff.empty()) continue;

        std::string key = "PI_TO_FF|" + port_name + "|" + ff;
        auto& e = edge_map[key];
        e.edge_type = "PI_TO_FF";
        e.port_name = port_name;
        e.ff_name = ff;
        e.ff_clock_net = io_ff_clock_cache.count(ff) ? io_ff_clock_cache[ff] : "";
        float slack = sta::delayAsFloat(pe->slack(sta_)) * SEC_TO_NS;
        if (!e.has_hold || slack < e.slack_hold_ns) {
          e.slack_hold_ns = slack;
          e.has_hold = true;
          pi_hold_ct++;
        }
      }
    }

    // --- REMOVAL (async hold equivalent, min) ---
    // Always extracted: removal violations appear in timing reports and CTS must
    // account for them via PI-hold-clip. If macro designs cause LP INFEASIBLE
    // due to massive removal edges, fix at LP solver level (not by hiding edges).
    {
      sta::PinSet* from_pins = new sta::PinSet(network_);
      from_pins->insert(port_pin);
      sta::ExceptionFrom* from = sta_->makeExceptionFrom(
          from_pins, nullptr, nullptr, sta::RiseFallBoth::riseFall());
      sta::PathEndSeq path_ends = sta_->findPathEnds(
          from, nullptr, nullptr, true, nullptr,
          sta::MinMaxAll::min(), 10000, 1, true, false,
          -sta::INF, sta::INF, false, nullptr,
          false, false, false, true, false, false);
      for (sta::PathEnd* pe : path_ends) {
        if (!pe) continue;
        sta::Vertex* end_vertex = pe->vertex(sta_);
        if (!end_vertex) continue;
        sta::Pin* end_pin = end_vertex->pin();
        std::string ff = getFFName(end_pin);
        if (ff.empty()) continue;
        std::string key = "PI_TO_FF|" + port_name + "|" + ff;
        auto& e = edge_map[key];
        e.edge_type = "PI_TO_FF"; e.port_name = port_name; e.ff_name = ff;
        e.ff_clock_net = io_ff_clock_cache.count(ff) ? io_ff_clock_cache[ff] : "";
        float slack = sta::delayAsFloat(pe->slack(sta_)) * SEC_TO_NS;
        if (!e.has_hold || slack < e.slack_hold_ns) {
          e.slack_hold_ns = slack; e.has_hold = true; pi_removal_ct++;
        }
      }
    }
    // --- RECOVERY (async setup equivalent, max) ---
    {
      sta::PinSet* from_pins = new sta::PinSet(network_);
      from_pins->insert(port_pin);
      sta::ExceptionFrom* from = sta_->makeExceptionFrom(
          from_pins, nullptr, nullptr, sta::RiseFallBoth::riseFall());
      sta::PathEndSeq path_ends = sta_->findPathEnds(
          from, nullptr, nullptr, true, nullptr,
          sta::MinMaxAll::max(), 10000, 1, true, false,
          -sta::INF, sta::INF, false, nullptr,
          false, false, true, false, false, false);
      for (sta::PathEnd* pe : path_ends) {
        if (!pe) continue;
        const sta::Path* path = pe->path();
        if (!path) continue;
        sta::Pin* end_pin = path->pin(sta_);
        std::string ff = getFFName(end_pin);
        if (ff.empty()) continue;
        std::string key = "PI_TO_FF|" + port_name + "|" + ff;
        auto& e = edge_map[key];
        e.edge_type = "PI_TO_FF"; e.port_name = port_name; e.ff_name = ff;
        e.ff_clock_net = io_ff_clock_cache.count(ff) ? io_ff_clock_cache[ff] : "";
        float slack = sta::delayAsFloat(pe->slack(sta_)) * SEC_TO_NS;
        if (!e.has_setup || slack < e.slack_setup_ns) {
          e.slack_setup_ns = slack; e.has_setup = true; pi_recovery_ct++;
        }
      }
    }
  }

  // ==================== FF→PO ====================
  for (auto* bterm : output_ports) {
    std::string port_name = bterm->getName();
    sta::Pin* port_pin = network_->dbToSta(bterm);
    if (!port_pin) continue;

    // --- SETUP (max) ---
    {
      sta::PinSet* to_pins = new sta::PinSet(network_);
      to_pins->insert(port_pin);
      sta::ExceptionTo* to = sta_->makeExceptionTo(
          to_pins, nullptr, nullptr,
          sta::RiseFallBoth::riseFall(), sta::RiseFallBoth::riseFall());

      sta::PathEndSeq path_ends = sta_->findPathEnds(
          nullptr, nullptr, to, true, nullptr,
          sta::MinMaxAll::max(), 10000, 1, true, false,
          -sta::INF, sta::INF, false, nullptr,
          true, false, false, false, false, false);

      for (sta::PathEnd* pe : path_ends) {
        if (!pe) continue;
        // For FF→PO, the startpoint is the FF
        const sta::Path* path = pe->path();
        if (!path) continue;
        // Walk to the start of the path to find the launch FF
        sta::PathExpanded expanded(path, sta_);
        if (expanded.size() == 0) continue;
        // path(0) = startpoint (launch FF), path(size()-1) = endpoint (PO)
        const sta::Path* start_path = expanded.path(0);
        if (!start_path) continue;
        sta::Pin* start_pin = start_path->pin(sta_);
        std::string ff = getFFName(start_pin);
        if (ff.empty()) continue;

        std::string key = "FF_TO_PO|" + port_name + "|" + ff;
        auto& e = edge_map[key];
        e.edge_type = "FF_TO_PO";
        e.port_name = port_name;
        e.ff_name = ff;
        e.ff_clock_net = io_ff_clock_cache.count(ff) ? io_ff_clock_cache[ff] : "";
        float slack = sta::delayAsFloat(pe->slack(sta_)) * SEC_TO_NS;
        if (!e.has_setup || slack < e.slack_setup_ns) {
          e.slack_setup_ns = slack;
          e.has_setup = true;
          po_setup_ct++;
        }
      }
    }

    // --- HOLD (min) ---
    {
      sta::PinSet* to_pins = new sta::PinSet(network_);
      to_pins->insert(port_pin);
      sta::ExceptionTo* to = sta_->makeExceptionTo(
          to_pins, nullptr, nullptr,
          sta::RiseFallBoth::riseFall(), sta::RiseFallBoth::riseFall());

      sta::PathEndSeq path_ends = sta_->findPathEnds(
          nullptr, nullptr, to, true, nullptr,
          sta::MinMaxAll::min(), 10000, 1, true, false,
          -sta::INF, sta::INF, false, nullptr,
          false, true, false, false, false, false);

      for (sta::PathEnd* pe : path_ends) {
        if (!pe) continue;
        const sta::Path* path = pe->path();
        if (!path) continue;
        sta::PathExpanded expanded(path, sta_);
        if (expanded.size() == 0) continue;
        // path(0) = startpoint (launch FF), path(size()-1) = endpoint (PO)
        const sta::Path* start_path = expanded.path(0);
        if (!start_path) continue;
        sta::Pin* start_pin = start_path->pin(sta_);
        std::string ff = getFFName(start_pin);
        if (ff.empty()) continue;

        std::string key = "FF_TO_PO|" + port_name + "|" + ff;
        auto& e = edge_map[key];
        e.edge_type = "FF_TO_PO";
        e.port_name = port_name;
        e.ff_name = ff;
        e.ff_clock_net = io_ff_clock_cache.count(ff) ? io_ff_clock_cache[ff] : "";
        float slack = sta::delayAsFloat(pe->slack(sta_)) * SEC_TO_NS;
        if (!e.has_hold || slack < e.slack_hold_ns) {
          e.slack_hold_ns = slack;
          e.has_hold = true;
          po_hold_ct++;
        }
      }
    }
  }

  // Convert map to vector
  result.reserve(edge_map.size());
  for (auto& [key, edge] : edge_map) {
    result.push_back(std::move(edge));
  }

  logger_->info(utl::CTS, 602,
                "3D-CTS IO: extracted {} IO edges "
                "(PI→FF: setup={}, hold={}, removal={}, recovery={}, "
                "FF→PO: setup={}, hold={})",
                result.size(), pi_setup_ct, pi_hold_ct,
                pi_removal_ct, pi_recovery_ct, po_setup_ct, po_hold_ct);

  return result;
}

void VerilogFFExtractor::writeIOCSV(const std::vector<IOEdgeVerilog>& edges,
                                    const std::string& output_file)
{
  std::ofstream file(output_file);
  if (!file.is_open()) {
    logger_->error(utl::CTS, 603, "Cannot open IO output file: {}", output_file);
  }

  // Same CSV format as Tcl extract_io_timing_edges_phase1a for cts_skew_lp.py
  file << "edge_type,port_name,ff_name,slack_setup_ns,slack_hold_ns,ff_clock_net\n";

  for (const auto& e : edges) {
    file << e.edge_type << "," << e.port_name << "," << e.ff_name << ",";
    if (e.has_setup) {
      file << e.slack_setup_ns;
    }
    file << ",";
    if (e.has_hold) {
      file << e.slack_hold_ns;
    }
    file << "," << e.ff_clock_net << "\n";
  }

  file.close();

  // Count by type
  int pi_ct = 0, po_ct = 0;
  for (const auto& e : edges) {
    if (e.edge_type == "PI_TO_FF") pi_ct++;
    else po_ct++;
  }

  logger_->info(utl::CTS, 604,
                "3D-CTS IO: wrote {} IO edges to {} (PI→FF: {}, FF→PO: {})",
                edges.size(), output_file, pi_ct, po_ct);
}

// ============================================================================
// BUF_MACRO: ODB-based register + macro sink collection.
// Phase 1: Iterates all dbInst, checks Liberty hasSequentials() for flip-flops.
// Also collects isBlock() macro instances (SRAMs).
//   SRAM blackboxes have Liberty timing arcs but no ff() groups, so
//   hasSequentials()=false. Physical CTS includes them via separateMacroRegSinks().
//   LP universe must match to enforce PI→macro and FF→macro hold constraints.
//   Without Phase 2, macro sinks get lp_target=0 and hold_budget=1.0ns (fail-open).
// ============================================================================
void VerilogFFExtractor::collectRegistersFromODB()
{
  if (!block_ || !network_) {
    logger_->warn(utl::CTS, 610,
                  "BUF_MACRO: block or network not available, cannot collect registers from ODB");
    return;
  }

  odb_registers_.clear();

  for (auto* inst : block_->getInsts()) {
    // Convert ODB inst to STA instance to access Liberty cell
    sta::Instance* sta_inst = network_->dbToSta(inst);
    if (!sta_inst) continue;

    sta::LibertyCell* lib_cell = network_->libertyCell(sta_inst);
    if (!lib_cell) continue;

    // hasSequentials() returns true for DFF, SDFF, DFFH, latch, etc.
    // Macros (SRAM, etc.) typically do NOT have hasSequentials() = true
    // because they are modeled as blackbox with timing arcs, not sequential cells.
    if (lib_cell->hasSequentials()) {
      std::string name = inst->getName();
      odb_registers_[name] = inst;
    }
  }

  const int seq_count = odb_registers_.size();

  // Also collect macro (block) instances that are
  // clock sinks. SRAM blackboxes have Liberty timing arcs but no ff() groups,
  // so hasSequentials() returns false. The physical CTS correctly includes them
  // via the same macro classifier used in separateMacroRegSinks(). 3D macro
  // views are not always consistent on dbInst::isBlock(), so we also fall back
  // to the master class. The LP universe must match
  // to enforce PI→macro and FF→macro hold constraints. Without this, macro sinks
  // get lp_target=0 and hold_budget=1.0ns (fail-open), causing severe hold
  // violations in macro designs (swerv_wrapper: hold 9.2x worse, ariane133: 3.3x).
  // Non-macro designs: both checks stay false for standard cells → no impact.
  int macro_count = 0;
  for (auto* inst : block_->getInsts()) {
    if (!isMacroBlockInst(inst)) continue;
    if (!inst->isPlaced()) continue;
    std::string name = inst->getName();
    if (odb_registers_.count(name)) continue;  // already collected as sequential
    odb_registers_[name] = inst;
    ++macro_count;
  }

  logger_->info(utl::CTS, 611,
                "BUF_MACRO: collected {} instances from ODB ({} sequential + {} macro block)",
                odb_registers_.size(), seq_count, macro_count);
}

// Per-clock-net LP: trace FF's clock input pin back through clock buffer chain
// to find the root clock net name (driven by a BTerm/port).
// This identifies which TritonCTS clock tree the FF belongs to.
// Example: FF clock pin -> clknet_leaf_5_clk -> clkbuf_3_clk -> clkbuf_0_clk -> clk (port)
// Returns "clk". For multi-clock designs (IBEX): some FFs trace to "clk",
// others to "clk_i_regs", etc.
std::string VerilogFFExtractor::traceClockNetRoot(odb::dbInst* ff_inst) const
{
  if (!ff_inst || !block_) return "";

  // Find the clock input ITerm (sigType == CLOCK)
  odb::dbITerm* clk_iterm = nullptr;
  for (odb::dbITerm* iterm : ff_inst->getITerms()) {
    if (iterm->getSigType() == odb::dbSigType::CLOCK) {
      clk_iterm = iterm;
      break;
    }
  }
  if (!clk_iterm) return "";

  odb::dbNet* net = clk_iterm->getNet();
  if (!net) return "";

  // Walk backward: net -> driver inst -> if single-output buffer, get input net -> repeat
  // Stop when net is driven by a BTerm (port) or we exceed max depth (safety).
  constexpr int MAX_DEPTH = 50;
  for (int depth = 0; depth < MAX_DEPTH; ++depth) {
    // Check if this net is driven by a port (BTerm)
    for (odb::dbBTerm* bterm : net->getBTerms()) {
      if (bterm->getIoType() == odb::dbIoType::INPUT) {
        return net->getName();  // Found root clock net
      }
    }

    // Find the driver instance of this net
    odb::dbITerm* driver = net->getFirstOutput();
    if (!driver) return net->getName();  // No driver, treat as root

    odb::dbInst* driver_inst = driver->getInst();
    if (!driver_inst) return net->getName();

    // Check if driver is a buffer/inverter (single input, single output)
    odb::dbITerm* input_iterm = nullptr;
    int input_count = 0;
    for (odb::dbITerm* iterm : driver_inst->getITerms()) {
      if (iterm->isInputSignal()) {
        input_iterm = iterm;
        input_count++;
      }
    }

    if (input_count != 1 || !input_iterm) {
      // Not a simple buffer (e.g., clock gater) — this is the effective root
      return net->getName();
    }

    // Continue tracing backward through the buffer's input net
    odb::dbNet* next_net = input_iterm->getNet();
    if (!next_net || next_net == net) return net->getName();  // Loop guard
    net = next_net;
  }

  return net->getName();  // Max depth reached
}

// ============================================================================
// BUF_MACRO: STA-based FF→FF edge extraction.
// For each register (from collectRegistersFromODB), calls STA findPathEnds
// to discover all FF→FF timing paths. This naturally traverses through:
//   - Combinational gates (standard cells)
//   - Macro blackboxes (Liberty timing arcs)
//   - Multi-level hierarchy
// Output format is identical to extractFFEdges() + fillTimingInfo().
//
// Key difference from Verilog-based:
//   - Verilog BFS: Q→net→gate→net→D (text parsing, misses macros)
//   - STA-based: findPathEnds from launch FF → all capture FFs (STA graph)
//
// Performance: O(N_ff × findPathEnds_cost). For AES ~534 FFs → 534 calls.
// Each call returns up to 10000 path ends (same as fillTimingInfo).
// ============================================================================
std::vector<FFEdgeVerilog> VerilogFFExtractor::extractFFEdgesViaSTA()
{
  std::vector<FFEdgeVerilog> edges;

  if (!sta_ || !network_ || !block_) {
    logger_->warn(utl::CTS, 612,
                  "BUF_MACRO: STA/network/block not available, cannot extract edges");
    return edges;
  }

  if (odb_registers_.empty()) {
    logger_->warn(utl::CTS, 613,
                  "BUF_MACRO: no registers found. Call collectRegistersFromODB() first.");
    return edges;
  }

  sta_->ensureGraph();
  sta_->ensureClkArrivals();
  sta_->updateTiming(false);

  constexpr float SEC_TO_NS = 1e9f;

  // Read EXTRACT_PATH_LIMIT env var — caps findPathEnds per launch FF.
  // No floor (unlike fillTimingInfoBatch): allows lower values to reduce LP scale.
  // Large macro designs (ariane133 25K FFs × 10K = 12M edges → LP INFEASIBLE).
  // Setting EXTRACT_PATH_LIMIT=500 → ~500K edges → LP feasible.
  // Timing fill functions (fillTimingInfo/Batch) use PATH_LIMIT with max(10000,...).
  int n_path_limit = 10000;  // backward-compat default
  const char* expl_env = std::getenv("EXTRACT_PATH_LIMIT");
  if (expl_env && std::atoi(expl_env) > 0) {
    n_path_limit = std::atoi(expl_env);
    logger_->info(utl::CTS, 605,
                  "extractFFEdgesViaSTA: EXTRACT_PATH_LIMIT={} per launch FF",
                  n_path_limit);
  }

  // Build ODB inst set for fast endpoint checking
  std::unordered_set<odb::dbInst*> register_set;
  for (const auto& [name, inst] : odb_registers_) {
    register_set.insert(inst);
  }

  // Helper: fill location + tier from ODB inst directly (no escapeForODB needed)
  auto fillFromODB = [&](odb::dbInst* inst, int& x, int& y, int& tier) {
    inst->getLocation(x, y);
    if (cts3dDb_) {
      tier = cts3dDb_->getInstTier(inst);
    } else {
      std::string master = inst->getMaster()->getName();
      tier = (master.find("_upper") != std::string::npos) ? 1 : 0;
    }
  };

  // Per-clock-net LP: build FF name -> root clock net cache.
  std::unordered_map<std::string, std::string> ff_clock_cache;
  for (const auto& [name, inst] : odb_registers_) {
    ff_clock_cache[name] = traceClockNetRoot(inst);
  }
  int n_traced = 0, n_untraceable = 0;
  for (const auto& [name, clk] : ff_clock_cache) {
    if (clk.empty()) n_untraceable++;
    else n_traced++;
  }
  logger_->info(utl::CTS, 616,
      "Per-clock-net LP: traced {} FFs to root clock nets ({} untraceable)",
      n_traced, n_untraceable);

  // Integer edge key: (src_idx << 32) | dst_idx — avoids string allocation/hash
  // for millions of edges. FF name → int index built from odb_registers_.
  std::unordered_map<std::string, int> ff_name_to_id;
  {
    int id = 0;
    for (const auto& [name, inst] : odb_registers_) {
      ff_name_to_id[name] = id++;
    }
  }
  // Also map capture FFs (endpoint) by odb::dbInst* → id for O(1) lookup
  std::unordered_map<odb::dbInst*, int> capture_odb_to_id;
  for (const auto& [name, inst] : odb_registers_) {
    capture_odb_to_id[inst] = ff_name_to_id[name];
  }

  // Deduplicate edges: integer key = (src_id << 32) | dst_id
  std::unordered_map<uint64_t, size_t> edge_map;

  int setup_found = 0;
  int hold_found = 0;
  int total_sta_calls = 0;
  int split_count = 0;

  // Build flat list of (name, odb_inst, sta_inst) for indexing
  struct SourceFF {
    std::string name;
    odb::dbInst* odb_inst;
    sta::Instance* sta_inst;
  };
  std::vector<SourceFF> all_sources;
  all_sources.reserve(odb_registers_.size());
  for (const auto& [name, inst] : odb_registers_) {
    sta::Instance* sta_inst = network_->dbToSta(inst);
    if (sta_inst) {
      all_sources.push_back({name, inst, sta_inst});
    }
  }

  // Process-level partitioning: EXTRACT_PARTITION="id/total" (e.g. "2/4")
  // Each partition processes a subset of source FFs.
  // Shell script launches N OpenROAD processes in parallel, then merges CSVs.
  int part_id = 0, part_total = 1;
  if (const char* p = std::getenv("EXTRACT_PARTITION")) {
    if (std::sscanf(p, "%d/%d", &part_id, &part_total) == 2
        && part_total > 1 && part_id >= 0 && part_id < part_total) {
      logger_->info(utl::CTS, 618,
          "extractFFEdgesViaSTA: partition {}/{} ({} of {} sources)",
          part_id, part_total,
          (static_cast<int>(all_sources.size()) + part_total - 1) / part_total,
          all_sources.size());
    } else {
      part_id = 0;
      part_total = 1;
    }
  }

  // Filter sources by partition
  std::vector<SourceFF> sources;
  if (part_total <= 1) {
    sources = std::move(all_sources);
  } else {
    for (int i = 0; i < static_cast<int>(all_sources.size()); ++i) {
      if (i % part_total == part_id) {
        sources.push_back(std::move(all_sources[i]));
      }
    }
  }

  // Pre-cache per-source ODB location + tier for O(1) lookup in batch results
  struct SourceInfo { int x; int y; int tier; };
  std::unordered_map<odb::dbInst*, SourceInfo> src_info_cache;
  for (const auto& s : sources) {
    SourceInfo si;
    fillFromODB(s.odb_inst, si.x, si.y, si.tier);
    src_info_cache[s.odb_inst] = si;
  }

  // O(1) launch FF lookup: odb::dbInst* -> source index
  std::unordered_map<odb::dbInst*, int> odb_to_src_idx;
  for (int i = 0; i < static_cast<int>(sources.size()); ++i) {
    odb_to_src_idx[sources[i].odb_inst] = i;
  }

  // Batch extraction disabled by default (K=1 = per-source).
  // Batch K>1 causes edge loss: dominant sources in batch monopolize path_limit,
  // dropping 69% edges vs per-source (swerv: 1.84M → 576K). Adaptive split
  // doesn't catch this because truncation occurs WITHIN the limit, not at it.
  // Forward scan + integer key still provide speedup with K=1.
  int batch_K = 1;
  if (const char* bk = std::getenv("EXTRACT_BATCH_SIZE")) {
    batch_K = std::max(1, std::atoi(bk));
  }

  logger_->info(utl::CTS, 617,
      "extractFFEdgesViaSTA: {} sources, batch_K={}, limit_per_src={}",
      sources.size(), batch_K, n_path_limit);

  // Helper: process findPathEnds results, extract edges, update edge_map.
  // Optimizations vs naive version:
  //   1. Launch FF recovery: forward scan (index 0→) instead of backward (size-1→0).
  //      Launch FF is at the START of PathExpanded → found in 1-3 hops, not 20-30.
  //   2. Integer edge key: uint64 instead of string concat+hash.
  //   3. Capture FF lookup via odb::dbInst* → int id (no getName()).
  auto processPathEnds = [&](const sta::PathEndSeq& path_ends,
                             const std::unordered_set<odb::dbInst*>& batch_set,
                             bool is_setup) {
    for (sta::PathEnd* path_end : path_ends) {
      if (!path_end) continue;

      const sta::Path* data_path = path_end->path();
      if (!data_path) continue;

      // Identify capture (endpoint) FF
      sta::Pin* end_pin = is_setup ? data_path->pin(sta_)
                                   : (path_end->vertex(sta_)
                                      ? path_end->vertex(sta_)->pin()
                                      : nullptr);
      if (!end_pin) continue;

      sta::Instance* end_inst = network_->instance(end_pin);
      if (!end_inst) continue;
      odb::dbInst* end_odb = network_->staToDb(end_inst);
      if (!end_odb || register_set.find(end_odb) == register_set.end()) continue;

      // Capture FF integer ID (O(1) lookup, no getName)
      auto cap_id_it = capture_odb_to_id.find(end_odb);
      if (cap_id_it == capture_odb_to_id.end()) continue;
      int dst_id = cap_id_it->second;

      // Identify launch FF: FORWARD scan (launch FF is near start of expanded path).
      // PathExpanded order: launch_clk → launch_FF → comb → capture_FF → capture_clk.
      // Forward scan finds launch in 1-3 hops vs 20-30 for backward scan.
      sta::PathExpanded expanded(data_path, sta_);
      odb::dbInst* launch_odb = nullptr;
      const int exp_sz = static_cast<int>(expanded.size());
      for (int pi = 0; pi < exp_sz; ++pi) {
        const sta::Path* pref = expanded.path(pi);
        if (!pref) continue;
        sta::Pin* pin = pref->pin(sta_);
        if (!pin) continue;
        sta::Instance* inst = network_->instance(pin);
        if (!inst) continue;
        odb::dbInst* odb_inst = network_->staToDb(inst);
        if (odb_inst && batch_set.count(odb_inst)) {
          launch_odb = odb_inst;
          break;
        }
      }
      if (!launch_odb) continue;

      // Source index + integer edge key
      auto src_it = odb_to_src_idx.find(launch_odb);
      if (src_it == odb_to_src_idx.end()) continue;
      int src_local_idx = src_it->second;
      const auto& src = sources[src_local_idx];

      auto src_id_it = ff_name_to_id.find(src.name);
      if (src_id_it == ff_name_to_id.end()) continue;
      int src_id = src_id_it->second;

      if (src_id == dst_id) continue;  // self-loop

      uint64_t key = (static_cast<uint64_t>(src_id) << 32)
                   | static_cast<uint64_t>(static_cast<uint32_t>(dst_id));

      float slack = sta::delayAsFloat(path_end->slack(sta_)) * SEC_TO_NS;
      float arrival = sta::delayAsFloat(data_path->arrival()) * SEC_TO_NS;
      float required = sta::delayAsFloat(path_end->requiredTime(sta_)) * SEC_TO_NS;

      auto it = edge_map.find(key);
      if (it == edge_map.end()) {
        // New edge — need names for CSV output
        FFEdgeVerilog edge;
        edge.from_ff = src.name;
        edge.to_ff = end_odb->getName();
        edge.from_clock_net = ff_clock_cache.count(src.name) ? ff_clock_cache[src.name] : "";
        edge.to_clock_net = ff_clock_cache.count(edge.to_ff) ? ff_clock_cache[edge.to_ff] : "";
        const auto& si = src_info_cache[launch_odb];
        edge.from_x = si.x; edge.from_y = si.y; edge.from_tier = si.tier;
        auto cap_si = src_info_cache.find(end_odb);
        if (cap_si != src_info_cache.end()) {
          edge.to_x = cap_si->second.x; edge.to_y = cap_si->second.y;
          edge.to_tier = cap_si->second.tier;
        } else {
          fillFromODB(end_odb, edge.to_x, edge.to_y, edge.to_tier);
        }
        if (is_setup) {
          edge.slack_max = slack; edge.arrival_max = arrival;
          edge.required_max = required; edge.has_timing = true;
          ++setup_found;
        } else {
          edge.slack_min = slack; edge.arrival_min = arrival;
          edge.required_min = required; edge.has_timing = true;
          ++hold_found;
        }
        edge_map[key] = edges.size();
        edges.push_back(edge);
      } else {
        auto& existing = edges[it->second];
        if (is_setup) {
          if (slack < existing.slack_max) {
            existing.slack_max = slack; existing.arrival_max = arrival;
            existing.required_max = required;
          }
        } else {
          if (existing.slack_min == 0.0f || slack < existing.slack_min) {
            existing.slack_min = slack; existing.arrival_min = arrival;
            existing.required_min = required;
          }
          ++hold_found;
        }
      }
    }
  };

  // Adaptive split batch extraction.
  // Calls findPathEnds with batch of K sources, limit=K*per_source_limit.
  // If result hits limit (truncation), splits batch in half and retries.
  // Base case K=1 is identical to original per-source extraction.
  std::function<void(int, int, bool)> batchExtract =
      [&](int start, int end_idx, bool is_setup) {
    int K = end_idx - start;
    if (K <= 0) return;

    int batch_limit = K * n_path_limit;

    // Build from-set for this batch
    sta::InstanceSet* from_insts = new sta::InstanceSet(network_);
    std::unordered_set<odb::dbInst*> batch_odb_set;
    for (int i = start; i < end_idx; ++i) {
      from_insts->insert(sources[i].sta_inst);
      batch_odb_set.insert(sources[i].odb_inst);
    }

    sta::ExceptionFrom* from = sta_->makeExceptionFrom(
        nullptr, nullptr, from_insts, sta::RiseFallBoth::riseFall());

    sta::PathEndSeq path_ends = sta_->findPathEnds(
        from, nullptr, nullptr, true, nullptr,
        is_setup ? sta::MinMaxAll::max() : sta::MinMaxAll::min(),
        batch_limit, 1, true, false,
        -sta::INF, sta::INF, false, nullptr,
        is_setup, !is_setup, false, false, false, false);

    ++total_sta_calls;

    int n_results = static_cast<int>(path_ends.size());

    if (n_results < batch_limit || K <= 1) {
      // No truncation (or single source) — process results
      processPathEnds(path_ends, batch_odb_set, is_setup);
    } else {
      // Truncation detected — split batch in half and retry
      ++split_count;
      int mid = start + K / 2;
      batchExtract(start, mid, is_setup);
      batchExtract(mid, end_idx, is_setup);
    }
  };

  // Run batched extraction for setup and hold
  int n_sources = static_cast<int>(sources.size());

  auto t0 = std::chrono::high_resolution_clock::now();

  for (int b_start = 0; b_start < n_sources; b_start += batch_K) {
    int b_end = std::min(b_start + batch_K, n_sources);
    batchExtract(b_start, b_end, true);   // setup
    batchExtract(b_start, b_end, false);  // hold

    // Progress log every 500 sources
    int processed = b_end;
    if (processed % 500 < batch_K || b_end == n_sources) {
      auto t1 = std::chrono::high_resolution_clock::now();
      double elapsed = std::chrono::duration<double>(t1 - t0).count();
      logger_->info(utl::CTS, 614,
          "extractFFEdgesViaSTA: {}/{} sources, {} edges, {} STA calls "
          "({} splits) in {:.1f}s",
          processed, n_sources, edges.size(), total_sta_calls,
          split_count, elapsed);
    }
  }

  auto t_end = std::chrono::high_resolution_clock::now();
  double total_elapsed = std::chrono::duration<double>(t_end - t0).count();

  logger_->info(utl::CTS, 615,
      "extractFFEdgesViaSTA: {} edges (setup={}, hold={}) from {} sources, "
      "{} STA calls ({} splits) in {:.1f}s",
      edges.size(), setup_found, hold_found, n_sources,
      total_sta_calls, split_count, total_elapsed);

  return edges;
}

}  // namespace cts
