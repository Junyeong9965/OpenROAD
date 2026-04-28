// SPDX-License-Identifier: BSD-3-Clause
// C++ LP-based buffer sizing solver
// Replaces Python buffer_sizing_lp.py with Liberty-based delay estimation.

#include "BufSizingLpSolver.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <regex>
#include <sstream>
#include <unordered_set>

#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "odb/db.h"
#include "sta/Corner.hh"
#include "sta/DcalcAnalysisPt.hh"
#include "sta/Liberty.hh"
#include "sta/MinMax.hh"
#include "sta/PathAnalysisPt.hh"
#include "sta/TableModel.hh"
#include "sta/PortDirection.hh"
#include "sta/TimingArc.hh"
#include "utl/Logger.h"

// OR-Tools GLOP
#include "ortools/linear_solver/linear_solver.h"

namespace cts {

BufSizingLpSolver::BufSizingLpSolver(odb::dbBlock* block,
                                       sta::dbSta* sta,
                                       utl::Logger* logger)
    : block_(block), sta_(sta), logger_(logger)
{
}

// Extract drive strength from cell name.
// Patterns: BUFx8_ASAP7... -> 8, BUF_X4_bottom -> 4, CLKBUF_X16... -> 16
int BufSizingLpSolver::getDriveStrength(const std::string& name)
{
  // Match BUFxN, BUF_XN, CLKBUFxN, CLKBUF_XN patterns (case-insensitive)
  std::regex rx(R"([xX](\d+))");
  std::smatch m;
  if (std::regex_search(name, m, rx)) {
    return std::stoi(m[1].str());
  }
  return 4;  // default
}

// Measure cell delay using Liberty NLDM table at given load capacitance.
// Same method as TechChar::charBufDelay (two-pass slew convergence).
double BufSizingLpSolver::measureCellDelay(sta::LibertyCell* cell,
                                            double load_cap)
{
  if (!cell) return 12.0;  // fallback 12ps

  sta::LibertyPort* in_port = nullptr;
  sta::LibertyPort* out_port = nullptr;
  cell->bufferPorts(in_port, out_port);
  if (!in_port || !out_port) return 12.0;

  sta::Corner* corner = sta_->cmdCorner();
  const sta::DcalcAnalysisPt* dcalc_ap
      = corner->findDcalcAnalysisPt(sta::MinMax::max());
  const sta::Pvt* pvt = dcalc_ap->operatingConditions();

  double max_delay = 0.0;
  for (sta::TimingArcSet* arcSet : cell->timingArcSets(in_port, out_port)) {
    for (sta::TimingArc* arc : arcSet->arcs()) {
      auto* model = dynamic_cast<sta::GateTimingModel*>(arc->model());
      if (!model) continue;
      sta::ArcDelay ad;
      sta::Slew sl;
      // Two-pass slew convergence (same as TechChar)
      model->gateDelay(pvt, 0.0, static_cast<float>(load_cap), false, ad, sl);
      model->gateDelay(pvt, sl, static_cast<float>(load_cap), false, ad, sl);
      double d = static_cast<double>(ad);
      if (d > max_delay) max_delay = d;
    }
  }

  return max_delay * 1e12;  // seconds -> ps
}

void BufSizingLpSolver::measureLibraryDelays()
{
  master_delays_.clear();
  auto* db_network = sta_->getDbNetwork();

  // Step 1: Enumerate ALL buffer masters in the design (from ODB)
  // This ensures we find all size variants (X1, X2, X4, X8, X16, etc.)
  std::unordered_set<std::string> master_names;
  for (const auto& buf : buffers_) {
    master_names.insert(buf.cell_name);
  }

  // Also find size variants by scanning ODB library for matching masters.
  // Include both BUF and INV cells — CTS uses INV as clkload cells.
  for (auto* lib : block_->getDb()->getLibs()) {
    for (auto* master : lib->getMasters()) {
      std::string mname = master->getName();
      std::string lower = mname;
      std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
      if (lower.find("buf") == std::string::npos
          && lower.find("inv") == std::string::npos) continue;
      master_names.insert(mname);
    }
  }

  // Step 2: Measure Liberty delay for each master.
  // Use bufferPorts() for buffers; for inverters, find input/output manually.
  for (const auto& name : master_names) {
    sta::LibertyCell* lib_cell = db_network->findLibertyCell(name.c_str());
    if (!lib_cell) continue;

    sta::LibertyPort* in_port = nullptr;
    sta::LibertyPort* out_port = nullptr;

    // Try buffer first
    lib_cell->bufferPorts(in_port, out_port);

    // If not a buffer, try to find single input/output (inverter)
    if (!in_port || !out_port) {
      sta::LibertyCellPortIterator port_iter(lib_cell);
      while (port_iter.hasNext()) {
        sta::LibertyPort* port = port_iter.next();
        if (port->direction() == sta::PortDirection::input()
            && !port->isClock()) {
          in_port = port;
        } else if (port->direction() == sta::PortDirection::output()) {
          out_port = port;
        }
      }
    }
    if (!in_port) continue;

    double input_cap = static_cast<double>(in_port->capacitance());
    double delay_ps = measureCellDelay(lib_cell, input_cap);
    master_delays_[name] = delay_ps;
  }

  logger_->info(utl::CTS, 890,
      "BufSizing: measured Liberty delays for {} masters", master_delays_.size());
}

void BufSizingLpSolver::extractBuffers()
{
  buffers_.clear();
  ff_to_buf_.clear();

  // Clock buffer/inverter name patterns (same as buffer_sizing.tcl + clkload)
  auto isClockBuf = [](const std::string& name) {
    return name.find("clkbuf") != std::string::npos
        || name.find("clkload") != std::string::npos
        || name.find("ghtap") != std::string::npos
        || name.find("skew_buf") != std::string::npos
        || name.find("cts") != std::string::npos;
  };

  // Iterate all instances
  for (auto* inst : block_->getInsts()) {
    std::string inst_name = inst->getName();
    if (!isClockBuf(inst_name)) continue;

    auto* master = inst->getMaster();
    if (!master) continue;
    std::string cell_name = master->getName();

    // Must be a buffer (1 input, 1 output signal pin)
    bool has_input = false, has_output = false;
    for (auto* mterm : master->getMTerms()) {
      if (mterm->getSigType() == odb::dbSigType::SIGNAL
          || mterm->getSigType() == odb::dbSigType::CLOCK) {
        if (mterm->getIoType() == odb::dbIoType::INPUT) has_input = true;
        if (mterm->getIoType() == odb::dbIoType::OUTPUT) has_output = true;
      }
    }
    if (!has_input || !has_output) continue;

    BufInfo buf;
    buf.inst_name = inst_name;
    buf.cell_name = cell_name;
    buf.buf_idx = static_cast<int>(buffers_.size());
    buf.current_delay_ps = 0.0;  // filled after measureLibraryDelays

    // Find output net -> driven FFs
    for (auto* iterm : inst->getITerms()) {
      if (iterm->getIoType() != odb::dbIoType::OUTPUT) continue;
      auto* net = iterm->getNet();
      if (!net) continue;
      for (auto* sink_iterm : net->getITerms()) {
        if (sink_iterm == iterm) continue;
        auto* sink_inst = sink_iterm->getInst();
        if (!sink_inst) continue;
        auto* sink_master = sink_inst->getMaster();
        if (!sink_master) continue;
        // Check if sink is sequential (FF)
        auto* db_net2 = sta_->getDbNetwork();
        sta::LibertyCell* lib_cell = db_net2->findLibertyCell(
            std::string(sink_master->getName()).c_str());
        if (lib_cell && lib_cell->isBuffer()) continue;  // skip downstream buffers
        if (lib_cell && lib_cell->hasSequentials()) {
          buf.driven_ffs.push_back(sink_inst->getName());
        }
      }
    }

    // Map FFs to buffer index
    for (const auto& ff : buf.driven_ffs) {
      ff_to_buf_[ff] = buf.buf_idx;
    }

    buffers_.push_back(std::move(buf));
  }

  logger_->info(utl::CTS, 891,
      "BufSizing: extracted {} clock buffers, {} FF mappings",
      buffers_.size(), ff_to_buf_.size());
}

void BufSizingLpSolver::extractTimingEdges(int /*max_paths*/)
{
  edges_.clear();

  // Read timing edges from CSV (Tcl extract_timing_edges_simple writes this).
  // CSV format: launch_ff,capture_ff,slack_setup_ps,slack_hold_ps,...
  const char* csv_env = std::getenv("BUF_SIZING_TIMING_CSV");
  if (!csv_env) {
    logger_->warn(utl::CTS, 902,
        "BufSizing: BUF_SIZING_TIMING_CSV not set, no timing edges");
    return;
  }

  std::ifstream f(csv_env);
  if (!f.is_open()) {
    logger_->warn(utl::CTS, 903,
        "BufSizing: cannot open timing CSV {}", csv_env);
    return;
  }

  std::string line;
  std::getline(f, line);  // skip header
  while (std::getline(f, line)) {
    std::stringstream ss(line);
    std::string token;
    std::vector<std::string> tokens;
    while (std::getline(ss, token, ',')) tokens.push_back(token);
    if (tokens.size() < 4) continue;

    TimingEdge edge;
    edge.launch_ff = tokens[0];
    edge.capture_ff = tokens[1];
    edge.slack_setup_ps = std::stod(tokens[2]);
    edge.slack_hold_ps = std::stod(tokens[3]);
    edges_.push_back(edge);
  }

  logger_->info(utl::CTS, 892,
      "BufSizing: loaded {} timing edges from {}", edges_.size(), csv_env);
}

std::vector<SizeOption> BufSizingLpSolver::getSizeOptions(
    const std::string& cell_name)
{
  std::vector<SizeOption> options;

  // Tier-based size option discovery (matches Python buffer_sizing_lp.py approach).
  // In heterogeneous 3D (asap7_nangate45_3D), CTS uses "BUF_X4_upper" (alias cell
  // with output pin Z) but actual ASAP7 size variants are "BUFx2_ASAP7_75t_R_upper"
  // etc. (output pin Y). Regex-based family matching FAILS because BUF_X{n}_upper
  // and BUFx{n}_ASAP7... are different naming families.
  //
  // Tier-based approach returns ALL buffer masters of the same tier.
  // Pin name mismatch (Z vs Y) is handled at swap time by the Tcl iterative
  // script (pin-remap swap), not filtered here — the LP needs to see all
  // sizing candidates to make globally optimal decisions.
  //
  // Algorithm:
  //   1. Detect tier suffix (_upper or _bottom)
  //   2. Return all BUF/INV masters in master_delays_ with matching tier suffix
  //   3. Filter out fractional-strength cells (INVxp*) unless current is also xp

  std::string tier_suffix;
  if (cell_name.size() >= 6
      && cell_name.substr(cell_name.size() - 6) == "_upper") {
    tier_suffix = "_upper";
  } else if (cell_name.size() >= 7
             && cell_name.substr(cell_name.size() - 7) == "_bottom") {
    tier_suffix = "_bottom";
  } else {
    return options;  // no tier suffix → cannot determine size variants
  }

  // Pre-compute current cell properties
  std::string cur_lower = cell_name;
  std::transform(cur_lower.begin(), cur_lower.end(), cur_lower.begin(),
                 ::tolower);
  bool current_is_inv = (cur_lower.find("inv") != std::string::npos);
  bool current_is_frac = (cur_lower.find("xp") != std::string::npos);

  // Collect all buffer/inverter masters with matching tier suffix
  for (const auto& [master_name, delay_ps] : master_delays_) {
    // Must end with same tier suffix
    if (master_name.size() < tier_suffix.size()) continue;
    if (master_name.substr(master_name.size() - tier_suffix.size()) != tier_suffix)
      continue;

    // Must be a buffer-type cell (BUF or INV in name, case-insensitive)
    std::string lower = master_name;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    bool is_buf = (lower.find("buf") != std::string::npos);
    bool is_inv = (lower.find("inv") != std::string::npos);
    if (!is_buf && !is_inv) continue;

    // Don't mix BUF and INV families: buffer→buffer, inverter→inverter.
    // Swapping BUF to INV inverts the signal polarity → functionally wrong.
    if (current_is_inv && !is_inv) continue;
    if (!current_is_inv && is_inv) continue;

    // Skip fractional-strength cells (INVxp33, INVxp67) unless current is also xp
    bool candidate_is_frac = (lower.find("xp") != std::string::npos);
    if (candidate_is_frac && !current_is_frac) continue;

    options.push_back({master_name, delay_ps});
  }

  // Sort by delay descending (smallest buf = highest delay first)
  std::sort(options.begin(), options.end(),
            [](const SizeOption& a, const SizeOption& b) {
              return a.delay_ps > b.delay_ps;
            });

  return options;
}

void BufSizingLpSolver::loadSkewTargets(const std::string& csv_path)
{
  buf_skew_targets_.clear();
  std::ifstream f(csv_path);
  if (!f.is_open()) return;

  std::string line;
  std::getline(f, line);  // header
  std::unordered_map<std::string, double> ff_targets;

  while (std::getline(f, line)) {
    std::stringstream ss(line);
    std::string ff_name, target_str;
    std::getline(ss, ff_name, ',');
    std::getline(ss, target_str, ',');
    double target_ns = std::stod(target_str);
    if (std::abs(target_ns) > 1e-9) {
      ff_targets[ff_name] = target_ns * 1e3;  // ns -> ps
    }
  }

  // Aggregate to per-buffer targets (mean of driven FF targets)
  for (int k = 0; k < static_cast<int>(buffers_.size()); ++k) {
    double sum = 0;
    int count = 0;
    for (const auto& ff : buffers_[k].driven_ffs) {
      auto it = ff_targets.find(ff);
      if (it != ff_targets.end()) {
        sum += it->second;
        count++;
      }
    }
    if (count > 0) {
      buf_skew_targets_[k] = sum / count;
    }
  }

  logger_->info(utl::CTS, 893,
      "BufSizing: loaded {} per-FF targets, aggregated to {} buffer targets",
      ff_targets.size(), buf_skew_targets_.size());
}

double BufSizingLpSolver::getMasterDelay(const std::string& master_name) const
{
  auto it = master_delays_.find(master_name);
  return (it != master_delays_.end()) ? it->second : 12.0;
}

std::vector<SizingChange> BufSizingLpSolver::solve()
{
  results_.clear();
  const int N = static_cast<int>(buffers_.size());
  if (N == 0) {
    logger_->info(utl::CTS, 894, "BufSizing: no buffers to optimize");
    return results_;
  }

  // Set current delays from Liberty measurements
  for (auto& buf : buffers_) {
    auto it = master_delays_.find(buf.cell_name);
    buf.current_delay_ps = (it != master_delays_.end()) ? it->second : 12.0;
  }

  // Candidate list generation mode (replaces LP optimization).
  // Generate ALL swap candidates for each buffer. The Tcl iterative script
  // will try each one with STA check (accept/reject based on TNS delta).
  //
  // Rationale: LP intrinsic delay optimization doesn't capture wire delay
  // effects from output cap changes (e.g., BUF_X4_upper→BUFx2_ASAP7 swap
  // reduces output cap → wire delay change → TNS improvement). Only STA
  // greedy can capture this. LP's role is just candidate generation.
  //
  // Strategy: for each buffer with skew target > 0, try ALL size variants
  // (sorted by drive strength ascending = downsize first, matching Python).
  // Buffers without skew target are also included (STA may find benefit).

  // Build set of buffers with skew targets (prioritize these)
  std::unordered_set<int> has_target;
  for (const auto& [k, target] : buf_skew_targets_) {
    if (k < N && target > 0.1) {
      has_target.insert(k);
    }
  }

  // Generate candidates: target buffers first, then non-target
  int n_target_cand = 0, n_other_cand = 0;
  auto gen_candidates = [&](int k) {
    auto options = getSizeOptions(buffers_[k].cell_name);
    for (const auto& opt : options) {
      if (opt.cell_name == buffers_[k].cell_name) continue;  // skip same cell
      int old_s = getDriveStrength(buffers_[k].cell_name);
      int new_s = getDriveStrength(opt.cell_name);
      results_.push_back({
          buffers_[k].inst_name,
          buffers_[k].cell_name,
          opt.cell_name,
          new_s > old_s ? "upsize" : "downsize"
      });
    }
    return static_cast<int>(options.size()) - 1;  // excluding self
  };

  // Priority 1: buffers with skew targets (downsize likely beneficial)
  for (int k = 0; k < N; ++k) {
    if (has_target.count(k)) {
      n_target_cand += gen_candidates(k);
    }
  }

  // Priority 2: remaining buffers (STA may find benefit)
  for (int k = 0; k < N; ++k) {
    if (!has_target.count(k)) {
      n_other_cand += gen_candidates(k);
    }
  }

  int up = 0, down = 0;
  for (const auto& r : results_) {
    if (r.direction == "upsize") up++;
    else down++;
  }
  logger_->info(utl::CTS, 901,
      "BufSizing candidates: {} total ({} from target bufs, {} other, {} upsize, {} downsize)",
      results_.size(), n_target_cand, n_other_cand, up, down);

  return results_;
}

void BufSizingLpSolver::writeResultsCSV(const std::string& csv_path) const
{
  std::ofstream f(csv_path);
  f << "inst_name,old_cell,new_cell,direction\n";
  for (const auto& r : results_) {
    f << r.inst_name << "," << r.old_cell << ","
      << r.new_cell << "," << r.direction << "\n";
  }
}

}  // namespace cts
