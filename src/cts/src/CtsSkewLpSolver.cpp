// SPDX-License-Identifier: BSD-3-Clause
// C++ LP solver using OR-Tools GLOP
// Ported from Python cts_skew_lp.py (V43/V44/V51 formulation)

#include "CtsSkewLpSolver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>

#include "Cts3DDatabase.h"
#include "odb/db.h"
#include "ortools/linear_solver/linear_solver.h"
#include "sta/Liberty.hh"
#include "sta/Network.hh"
#include "sta/Sta.hh"
#include "utl/Logger.h"

namespace cts {

using operations_research::MPConstraint;
using operations_research::MPObjective;
using operations_research::MPSolver;
using operations_research::MPVariable;

// Log message IDs for CTS module (CTS-600 series for LP solver)
constexpr int CTS_LP_START = 600;
constexpr int CTS_LP_EDGES = 601;
constexpr int CTS_LP_VARS = 602;
constexpr int CTS_LP_RESULT = 603;
constexpr int CTS_LP_TARGETS = 604;
constexpr int CTS_LP_IO = 605;

CtsSkewLpSolver::CtsSkewLpSolver(odb::dbBlock* block,
                                 sta::dbSta* sta,
                                 sta::dbNetwork* network,
                                 utl::Logger* logger,
                                 Cts3DDatabase* cts3dDb)
    : block_(block),
      sta_(sta),
      network_(network),
      logger_(logger),
      cts3dDb_(cts3dDb)
{
}

bool CtsSkewLpSolver::isGhostEdge(float slack_max, float slack_min,
                                   float arr_max, float arr_min) const
{
  // Ghost edges have all-zero timing fields (STA couldn't find valid path).
  // These create contradictory hold constraints → LP INFEASIBLE.
  // Filter: all four values must be zero (avoids false positives).
  return std::abs(slack_max) < 1e-12f
      && std::abs(slack_min) < 1e-12f
      && std::abs(arr_max)   < 1e-12f
      && std::abs(arr_min)   < 1e-12f;
}

// Read propagated-clock CSV from Phase 1 sub-OpenROAD.
// Previous code used live STA (main session, ideal-clock) → only 24 violations vs 2189.
// Phase 1 sub-process runs balanced CTS + set_propagated_clock + extraction → CSV has
// correct propagated-clock slacks that reflect H-tree latency differences.
void CtsSkewLpSolver::readTimingGraphCSV(const std::string& csv_path,
                                          const std::string& clock_net_filter)
{
  logger_->info(utl::CTS, CTS_LP_START,
                "3D-CTS LP: Reading timing graph from CSV: {}", csv_path);

  std::ifstream infile(csv_path);
  if (!infile.is_open()) {
    logger_->error(utl::CTS, CTS_LP_START + 1,
                   "3D-CTS LP: Cannot open timing graph CSV: {}", csv_path);
    return;
  }

  // Parse header
  // Format: from_ff,to_ff,slack_max_ns,slack_min_ns,arrival_max_ns,arrival_min_ns,
  //         required_max_ns,required_min_ns,from_x,from_y,from_tier,to_x,to_y,to_tier
  std::string header;
  std::getline(infile, header);

  // Build FF set and edges from CSV
  struct RawEdge {
    std::string from_ff, to_ff;
    float slack_max, slack_min, arrival_max, arrival_min;
    int from_tier, to_tier;
    std::string from_clock_net, to_clock_net;  // Per-clock-net LP
  };
  std::vector<RawEdge> raw_edges;
  std::unordered_set<std::string> ff_set;

  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty()) continue;
    std::istringstream ss(line);
    std::string token;
    std::vector<std::string> fields;
    while (std::getline(ss, token, ',')) {
      fields.push_back(token);
    }
    if (fields.size() < 14) continue;

    // Handle "inf"/"-inf"/empty/NaN in timing graph CSV.
    // Python float("inf") = inf → edge kept (inf != 0.0 → not ghost → FF in set,
    // inf slack → eff_setup=inf → never violated → constraint harmless).
    // C++ must match: return actual infinity for "inf" strings so FF stays in set
    // and edge passes ghost filter. Inf-slack edges are later skipped when building
    // edges_ (no real timing constraint). Empty/NaN → 0.0 (ghost filter handles).
    // Previous approach (skip inf edge entirely) lost ~1100 FFs in nangate45_3D/ibex
    // → hold:0 → LP distortion → CKMeans topology destroyed → -188k ps regression.
    auto safe_stof = [](const std::string& s) -> float {
      if (s.empty()) return 0.0f;
      if (s == "inf" || s == "INF") {
        return std::numeric_limits<float>::infinity();
      }
      if (s == "-inf" || s == "-INF") {
        return -std::numeric_limits<float>::infinity();
      }
      if (s == "nan" || s == "NAN") {
        return 0.0f;
      }
      try {
        return std::stof(s);
      } catch (...) {
        return 0.0f;
      }
    };

    RawEdge e;
    e.from_ff = fields[0];
    e.to_ff = fields[1];
    e.slack_max = safe_stof(fields[2]);
    e.slack_min = safe_stof(fields[3]);
    e.arrival_max = safe_stof(fields[4]);
    e.arrival_min = safe_stof(fields[5]);
    e.from_tier = std::stoi(fields[10]);
    e.to_tier = std::stoi(fields[13]);

    // Per-clock-net LP: columns 14,15 (from_clock_net, to_clock_net)
    // Backward compatible: if CSV has only 14 columns, clock_net stays empty.
    if (fields.size() >= 16) {
      e.from_clock_net = fields[14];
      e.to_clock_net = fields[15];
    }

    // Keep ALL edges in raw_edges (including inf-slack).
    // Inf edges contribute FFs to ff_set but are skipped when building edges_.
    raw_edges.push_back(e);
  }

  // Build FF set AFTER ghost edge filtering.
  // Previously ff_set was populated before ghost filter → ghost-only FFs
  // (e.g., _16656__bottom) got included with t_max=100ps and a variable
  // in the LP, causing 535 FFs vs Python's 534. Now matches Python behavior.
  edges_.clear();
  int n_ghost = 0;
  int n_cross_clock = 0;
  for (const auto& e : raw_edges) {
    if (isGhostEdge(e.slack_max, e.slack_min, e.arrival_max, e.arrival_min)) {
      ++n_ghost;
      continue;
    }
    // Per-clock-net LP: skip edges not matching the clock net filter.
    // Both endpoints must belong to the target clock net.
    // Cross-clock edges (from_clock_net != to_clock_net) are always excluded when filtering.
    if (!clock_net_filter.empty()) {
      if (e.from_clock_net != clock_net_filter || e.to_clock_net != clock_net_filter) {
        ++n_cross_clock;
        continue;
      }
    }
    ff_set.insert(e.from_ff);
    ff_set.insert(e.to_ff);
  }

  // Build FF list and index map
  ff_list_.assign(ff_set.begin(), ff_set.end());
  std::sort(ff_list_.begin(), ff_list_.end());
  ff_to_idx_.clear();
  for (int i = 0; i < static_cast<int>(ff_list_.size()); ++i) {
    ff_to_idx_[ff_list_[i]] = i;
  }

  // Fill tiers from non-ghost edges
  ff_tiers_.resize(ff_list_.size(), 0);
  for (const auto& e : raw_edges) {
    if (isGhostEdge(e.slack_max, e.slack_min, e.arrival_max, e.arrival_min)) {
      continue;
    }
    auto it = ff_to_idx_.find(e.from_ff);
    if (it != ff_to_idx_.end()) ff_tiers_[it->second] = e.from_tier;
    it = ff_to_idx_.find(e.to_ff);
    if (it != ff_to_idx_.end()) ff_tiers_[it->second] = e.to_tier;
  }

  // Convert to LpTimingEdge
  // Skip ghost edges and inf-slack edges (unconstrained paths).
  // Inf edges already contributed FFs to ff_set above; they carry no
  // real timing constraint, so they must not become LP constraints.
  for (const auto& e : raw_edges) {
    if (isGhostEdge(e.slack_max, e.slack_min, e.arrival_max, e.arrival_min)) {
      continue;
    }
    if (std::isinf(e.slack_max) || std::isinf(e.slack_min)) {
      continue;
    }
    // Per-clock-net LP: skip cross-clock and non-matching edges
    if (!clock_net_filter.empty()) {
      if (e.from_clock_net != clock_net_filter || e.to_clock_net != clock_net_filter) {
        continue;
      }
    }
    auto it_i = ff_to_idx_.find(e.from_ff);
    auto it_j = ff_to_idx_.find(e.to_ff);
    if (it_i == ff_to_idx_.end() || it_j == ff_to_idx_.end()) {
      continue;
    }
    LpTimingEdge le;
    le.from_idx = it_i->second;
    le.to_idx = it_j->second;
    le.slack_setup_ns = e.slack_max;
    le.slack_hold_ns = e.slack_min;
    le.cross_tier = (e.from_tier != e.to_tier);
    // Hold-ghost: STA found no hold (min-delay) path for this edge.
    // slack_min=0 AND arrival_min=0 means unconstrained, not "hold slack = 0".
    le.hold_ghost = (std::abs(e.slack_min) < 1e-12f && std::abs(e.arrival_min) < 1e-12f);
    edges_.push_back(le);
  }

  logger_->info(utl::CTS, CTS_LP_EDGES,
                "3D-CTS LP: {} FFs, {} edges ({} ghost, {} cross-clock filtered, {} cross-tier){}",
                ff_list_.size(), edges_.size(), n_ghost, n_cross_clock,
                std::count_if(edges_.begin(), edges_.end(),
                              [](const LpTimingEdge& e) { return e.cross_tier; }),
                clock_net_filter.empty() ? "" : " [clock: " + clock_net_filter + "]");
}

// Read IO edges from Phase 1 CSV (propagated-clock).
// Format: edge_type,port_name,ff_name,slack_setup_ns,slack_hold_ns
void CtsSkewLpSolver::readIOEdgesCSV(const std::string& io_csv_path,
                                      const std::string& clock_net_filter)
{
  io_edges_.clear();

  if (io_csv_path.empty()) {
    logger_->info(utl::CTS, CTS_LP_IO,
                  "3D-CTS LP: IO edges: 0 (no IO CSV path)");
    return;
  }

  std::ifstream infile(io_csv_path);
  if (!infile.is_open()) {
    logger_->warn(utl::CTS, CTS_LP_IO,
                  "3D-CTS LP: Cannot open IO CSV: {} — skipping IO edges", io_csv_path);
    return;
  }

  std::string header;
  std::getline(infile, header);

  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty()) continue;
    std::istringstream ss(line);
    std::string token;
    std::vector<std::string> fields;
    while (std::getline(ss, token, ',')) {
      fields.push_back(token);
    }
    if (fields.size() < 5) continue;

    const std::string& edge_type = fields[0];
    const std::string& ff_name = fields[2];

    // Per-clock-net LP: parse ff_clock_net (column 5, optional)
    std::string ff_clock_net;
    if (fields.size() >= 6) {
      ff_clock_net = fields[5];
    }

    // Skip IO edges for FFs not on the target clock net
    if (!clock_net_filter.empty() && ff_clock_net != clock_net_filter) {
      continue;
    }

    // Handle empty/inf/nan fields in IO edge CSV.
    // Tcl extraction writes empty strings when STA find_timing_paths doesn't
    // find a valid path for certain (port, FF) pairs. std::stof("") throws
    // → LP solver crashes → all targets=0 → balanced CTS only.
    // nangate45_3D/ibex: Setup TNS -210,266ps (6.8x worse than Pin3D) was
    // caused by this crash. Apply same safe_stof pattern as timing graph parser.
    auto safe_stof_io = [](const std::string& s) -> std::pair<bool, float> {
      if (s.empty()) return {false, 0.0f};
      if (s == "inf" || s == "INF" || s == "-inf" || s == "-INF" ||
          s == "nan" || s == "NAN") {
        return {false, 0.0f};
      }
      try {
        return {true, std::stof(s)};
      } catch (...) {
        return {false, 0.0f};
      }
    };

    auto [valid_setup, slack_setup] = safe_stof_io(fields[3]);
    auto [valid_hold, slack_hold] = safe_stof_io(fields[4]);
    if (!valid_setup || !valid_hold) continue;  // Skip incomplete IO edges

    auto it = ff_to_idx_.find(ff_name);
    if (it == ff_to_idx_.end()) continue;

    if (edge_type == "PI_TO_FF") {
      // PI hold (most important — prevents over-delivery on input paths)
      LpIOEdge le;
      le.ff_idx = it->second;
      le.slack_ns = slack_hold;
      le.type = LpIOEdge::PI_HOLD;
      io_edges_.push_back(le);
      // PI setup
      LpIOEdge le2;
      le2.ff_idx = it->second;
      le2.slack_ns = slack_setup;
      le2.type = LpIOEdge::PI_SETUP;
      io_edges_.push_back(le2);
    } else if (edge_type == "FF_TO_PO") {
      LpIOEdge le;
      le.ff_idx = it->second;
      le.slack_ns = slack_setup;
      le.type = LpIOEdge::PO_SETUP;
      io_edges_.push_back(le);
      LpIOEdge le2;
      le2.ff_idx = it->second;
      le2.slack_ns = slack_hold;
      le2.type = LpIOEdge::PO_HOLD;
      io_edges_.push_back(le2);
    }
  }

  int n_pi_hold = 0, n_pi_setup = 0, n_po_setup = 0, n_po_hold = 0;
  for (const auto& io : io_edges_) {
    switch (io.type) {
      case LpIOEdge::PI_HOLD:  ++n_pi_hold;  break;
      case LpIOEdge::PI_SETUP: ++n_pi_setup; break;
      case LpIOEdge::PO_SETUP: ++n_po_setup; break;
      case LpIOEdge::PO_HOLD:  ++n_po_hold;  break;
    }
  }
  logger_->info(utl::CTS, CTS_LP_IO,
                "3D-CTS LP: IO edges: {} total (pi_hold:{} pi_setup:{} po_setup:{} po_hold:{})",
                io_edges_.size(), n_pi_hold, n_pi_setup, n_po_setup, n_po_hold);
}

void CtsSkewLpSolver::computeBounds()
{
  // Removed hold-slack-based t_max clipping.
  // Previous code clipped t_max using min_hold_as_capture, which:
  //   (a) duplicates the LP's own hold constraints (double-counting)
  //   (b) over-restricts the solution space — tighter bounds cause the LP
  //       to miss setup improvement opportunities
  //   (c) differs from Python solver which uses --bounds-csv (physical
  //       achievability from ClockLatencyEstimator) or plain max_skew
  // The LP's hard hold constraints (a_j - a_i <= eff_hold - margin) already
  // enforce hold safety. PI-hold-clip (applyPiHoldClip) handles PI→FF hold.
  // Now simply initializes t_max_ with max_skew, matching Python behavior.

  const int n_ffs = static_cast<int>(ff_list_.size());
  t_max_.resize(n_ffs, params_.max_skew);

  // delivery_max clamp removed — LP uses max_skew only.
  // CTS builder caps at max_depth × d_buf at delivery time.

  logger_->info(utl::CTS, CTS_LP_IO + 2,
                "3D-CTS LP: Bounds: {} FFs, t_max={:.1f}ps (no delivery clamp)",
                n_ffs, params_.max_skew * 1e3);
}

void CtsSkewLpSolver::applyPiHoldClip()
{
  // PI-hold-clip: for FFs with PI→FF hold/removal path, reduce t_max
  // to prevent LP from assigning arrival that violates input-path hold.
  // This enforces: delivery_bound_per_FF ≤ min(global_delivery_max, PI_hold_slack).
  //
  // Covers both synchronous hold (data D pin) and asynchronous removal
  // (reset/set RESETN/SETN pin) paths — both are PI_HOLD type in io_edges_
  // after the V54 extraction fix that adds removal pass.
  //
  // For each PI→FF hold/removal edge with slack > 0:
  //   t_max[ff] = min(t_max[ff], slack - sigma_pi)
  // For each PI→FF hold/removal edge with slack <= 0:
  //   t_max[ff] = 0  (cannot add any delay)

  int n_zeroed = 0;
  int n_clipped = 0;
  int n_delivery_limited = 0;  // FFs where PI clip < delivery_max
  for (const auto& io : io_edges_) {
    if (io.type != LpIOEdge::PI_HOLD) continue;
    int idx = io.ff_idx;
    if (idx < 0 || idx >= static_cast<int>(t_max_.size())) continue;

    const float prev_tmax = t_max_[idx];
    float slack = io.slack_ns;
    if (slack <= 0) {
      t_max_[idx] = 0.0f;
      ++n_zeroed;
    } else {
      float clipped = std::max(0.0f, slack - params_.sigma_pi);
      if (clipped < t_max_[idx]) {
        t_max_[idx] = clipped;
        ++n_clipped;
      }
    }
    // Count FFs where PI-hold-clip reduced t_max below delivery_max
    // (i.e., PI hold slack is the binding constraint, not physical delivery)
    if (params_.delivery_max > 0 && t_max_[idx] < prev_tmax) {
      ++n_delivery_limited;
    }
  }

  if (n_zeroed > 0 || n_clipped > 0) {
    logger_->info(utl::CTS, CTS_LP_IO + 1,
                  "3D-CTS LP: PI-hold-clip: {} zeroed, {} clipped, "
                  "{} delivery-limited (sigma_pi={:.1f}ps)",
                  n_zeroed, n_clipped, n_delivery_limited,
                  params_.sigma_pi * 1e3);
  }
}

LpResult CtsSkewLpSolver::solveTns()
{
  const int n_ffs = static_cast<int>(ff_list_.size());
  if (n_ffs == 0) {
    result_.status = "no_ffs";
    return result_;
  }

  const float sigma = params_.sigma_local;

  // Compute effective slacks and classify edges
  struct SetupEdge { int i; int j; float eff_slack; };
  struct HoldEdge  { int i; int j; float eff_slack; };

  std::vector<SetupEdge> setup_edges;
  std::vector<HoldEdge>  hold_edges;

  int n_hold_ghost = 0;
  for (const auto& e : edges_) {
    float eff_setup = e.slack_setup_ns - 2.0f * sigma;  // sigma_i + sigma_j
    setup_edges.push_back({e.from_idx, e.to_idx, eff_setup});

    // Skip hold-ghost edges: STA found no hold (min-delay) path.
    // slack_min=0 AND arrival_min=0 means unconstrained, not "hold slack = 0".
    // Without this filter, sigma correction turns slack=0 → eff_hold=-10ps →
    // fake negative-rhs hold constraints → massive cycles → LP INFEASIBLE.
    // IBEX: 483K/1.14M edges (42%) are hold-ghost.
    if (e.hold_ghost) {
      ++n_hold_ghost;
      continue;
    }
    float eff_hold = e.slack_hold_ns - 2.0f * sigma;
    hold_edges.push_back({e.from_idx, e.to_idx, eff_hold});
  }

  int n_neg_setup = std::count_if(setup_edges.begin(), setup_edges.end(),
                                   [](const SetupEdge& e) { return e.eff_slack < 0; });

  logger_->info(utl::CTS, CTS_LP_VARS,
                "3D-CTS LP: {} setup edges ({} violating), {} hold edges",
                setup_edges.size(), n_neg_setup, hold_edges.size());

  if (n_neg_setup == 0) {
    // No violations — all targets = 0
    result_.status = "no_violations";
    for (const auto& ff : ff_list_) {
      result_.arrivals[ff] = 0.0;
    }
    result_.n_ffs = n_ffs;
    return result_;
  }

  // Prune setup edges: skip if eff_slack > t_max_i + t_max_j
  struct PrunedSetup { int i; int j; float eff_slack; };
  std::vector<PrunedSetup> pruned_setup;
  for (const auto& e : setup_edges) {
    if (e.eff_slack <= t_max_[e.i] + t_max_[e.j]) {
      pruned_setup.push_back({e.i, e.j, e.eff_slack});
    }
  }

  // Separate hard vs soft IO constraints
  std::vector<std::pair<int, float>> hard_io;   // (ff_idx, corrected_slack)
  std::vector<std::tuple<int, float, LpIOEdge::Type>> soft_io;  // (ff_idx, slack, type)

  for (const auto& io : io_edges_) {
    if (io.type == LpIOEdge::PI_HOLD && params_.hard_pi_hold) {
      float corrected = std::max(0.0f, io.slack_ns);
      if (corrected > 0) {
        hard_io.push_back({io.ff_idx, corrected});
      } else {
        soft_io.push_back({io.ff_idx, io.slack_ns, io.type});
      }
    } else {
      soft_io.push_back({io.ff_idx, io.slack_ns, io.type});
    }
  }
  const int n_io_soft = static_cast<int>(soft_io.size());

  // ===== CG mode check =====
  // Constraint Generation: setup edges added only when violated.
  // Effective for small/medium designs (AES 534 FFs: 15% edges added = 85% LP reduction).
  // Ineffective for large designs (IBEX 10K FFs: 99% edges added = no pruning benefit).
  // Auto-detect: CG ON if n_ffs <= 5000, OFF otherwise.
  // CTS_ENABLE_CG env var overrides auto-detect (0=force OFF, 1=force ON).
  constexpr int CG_FF_THRESHOLD = 5000;
  bool enable_cg = (n_ffs <= CG_FF_THRESHOLD);
  if (const char* e = std::getenv("CTS_ENABLE_CG")) {
    enable_cg = (std::atoi(e) != 0);
    logger_->info(utl::CTS, 703,
        "3D-CTS LP: CTS_ENABLE_CG={} (env override, n_ffs={})",
        enable_cg ? 1 : 0, n_ffs);
  } else {
    logger_->info(utl::CTS, 705,
        "3D-CTS LP: CG auto-detect n_ffs={} {} threshold={} → CG {}",
        n_ffs, n_ffs <= CG_FF_THRESHOLD ? "<=" : ">",
        CG_FF_THRESHOLD, enable_cg ? "ON" : "OFF");
  }

  // TopK: only when CG is OFF (legacy path)
  const int topk = enable_cg ? 0 : params_.endpoint_topk;
  if (topk > 0 && !pruned_setup.empty()) {
    const int pre_topk = static_cast<int>(pruned_setup.size());
    std::unordered_map<int, std::vector<int>> ep_indices;
    for (int idx = 0; idx < static_cast<int>(pruned_setup.size()); ++idx) {
      ep_indices[pruned_setup[idx].j].push_back(idx);
    }
    std::vector<bool> keep(pruned_setup.size(), false);
    for (auto& [ep_j, indices] : ep_indices) {
      if (static_cast<int>(indices.size()) <= topk) {
        for (int idx : indices) keep[idx] = true;
      } else {
        std::partial_sort(indices.begin(), indices.begin() + topk, indices.end(),
            [&](int a, int b) {
              return pruned_setup[a].eff_slack < pruned_setup[b].eff_slack;
            });
        for (int k = 0; k < topk; ++k) keep[indices[k]] = true;
      }
    }
    std::vector<PrunedSetup> topk_setup;
    topk_setup.reserve(pruned_setup.size());
    for (int idx = 0; idx < static_cast<int>(pruned_setup.size()); ++idx) {
      if (keep[idx]) topk_setup.push_back(pruned_setup[idx]);
    }
    const int post_topk = static_cast<int>(topk_setup.size());
    if (post_topk < pre_topk) {
      logger_->info(utl::CTS, 690,
          "3D-CTS LP: Endpoint-TopK={}: setup edges {} -> {} ({} endpoints)",
          topk, pre_topk, post_topk,
          static_cast<int>(ep_indices.size()));
    }
    pruned_setup = std::move(topk_setup);
  }

  // ===== Create OR-Tools GLOP solver =====
  std::unique_ptr<MPSolver> solver(
      MPSolver::CreateSolver("GLOP"));
  if (!solver) {
    logger_->error(utl::CTS, CTS_LP_RESULT,
                   "3D-CTS LP: Failed to create GLOP solver");
    result_.status = "solver_creation_failed";
    return result_;
  }

  const double inf = solver->infinity();

  // a_i: per-FF arrival in [0, t_max_i]
  std::vector<MPVariable*> a_vars(n_ffs);
  for (int i = 0; i < n_ffs; ++i) {
    a_vars[i] = solver->MakeNumVar(0.0, t_max_[i],
                                    "a_" + std::to_string(i));
  }

  // w_m: IO soft penalty variables, >= 0
  std::vector<MPVariable*> w_vars(n_io_soft);
  for (int m = 0; m < n_io_soft; ++m) {
    w_vars[m] = solver->MakeNumVar(0.0, inf, "w_" + std::to_string(m));
  }

  // W: WNS minimax variable
  MPVariable* W_var = nullptr;
  if (params_.gamma_wns > 0 || params_.enable_2phase) {
    W_var = solver->MakeNumVar(0.0, inf, "W");
  }

  // ===== Setup edge variables (compact formulation) =====
  // Compact: no per-edge v_k. V_j absorbs max violation directly.
  //   a_i - a_j - V_j ≤ eff_slack  (equivalent to v_k + linkage, 95% fewer vars)
  std::unordered_map<int, MPVariable*> Vj_map;  // j (FF idx) -> Vj var

  int n_setup_ct = 0;
  int n_wns_ct = 0;

  if (!enable_cg) {
    // Non-CG: create all V_j and setup constraints upfront
    for (const auto& ps : pruned_setup) {
      if (Vj_map.find(ps.j) == Vj_map.end()) {
        auto* Vj = solver->MakeNumVar(0.0, inf, "Vj_" + std::to_string(ps.j));
        Vj_map[ps.j] = Vj;
      }
    }
    // Compact setup: a_i - a_j - V_j ≤ eff_slack
    for (const auto& ps : pruned_setup) {
      MPConstraint* ct = solver->MakeRowConstraint(-inf, ps.eff_slack, "");
      ct->SetCoefficient(a_vars[ps.i], 1.0);
      ct->SetCoefficient(a_vars[ps.j], -1.0);
      ct->SetCoefficient(Vj_map.at(ps.j), -1.0);
      ++n_setup_ct;
    }
    // WNS minimax: V_j - W ≤ 0
    if (W_var) {
      for (auto& [j, Vj] : Vj_map) {
        MPConstraint* ct = solver->MakeRowConstraint(-inf, 0.0, "");
        ct->SetCoefficient(Vj, 1.0);
        ct->SetCoefficient(W_var, -1.0);
        ++n_wns_ct;
      }
    }
  }
  // CG mode: Vj_map starts empty — populated in CG loop below.

  // ===== Objective =====
  MPObjective* objective = solver->MutableObjective();
  objective->SetMinimization();

  // Phase B objective setter (compact: no v_k, only V_j)
  auto setPhaseBObjective = [&]() {
    for (int i = 0; i < n_ffs; ++i)
      objective->SetCoefficient(a_vars[i], params_.lambda_reg);
    for (auto& [j, Vj] : Vj_map)
      objective->SetCoefficient(Vj, 1.0);
    for (int m = 0; m < n_io_soft; ++m)
      objective->SetCoefficient(w_vars[m], params_.weight_io);
    if (W_var) objective->SetCoefficient(W_var, params_.gamma_wns);
  };

  auto setPhaseAObjective = [&]() {
    if (W_var) objective->SetCoefficient(W_var, 1.0);
    for (int i = 0; i < n_ffs; ++i)
      objective->SetCoefficient(a_vars[i], params_.lambda_reg * 1e-4f);
    for (auto& [j, Vj] : Vj_map)
      objective->SetCoefficient(Vj, 0.0);
    for (int m = 0; m < n_io_soft; ++m)
      objective->SetCoefficient(w_vars[m], 0.0);
  };

  // ===== Non-setup constraints (built upfront for both CG and non-CG) =====
  int n_hold_ct = 0;
  int n_hard_ct = 0;
  int n_io_ct = 0;

  // Hold: a_j - a_i ≤ eff_hold - hold_margin  (HARD)
  // Step 1: Single-edge skip (redundant + structurally infeasible)
  struct KeptHold { int i; int j; float rhs; };
  std::vector<KeptHold> hold_kept;
  int n_hold_skip_infeasible = 0;
  for (const auto& he : hold_edges) {
    float rhs = he.eff_slack - params_.hold_margin;
    if (rhs > t_max_[he.i] + t_max_[he.j]) continue;
    if (rhs < 0.0f && (t_max_[he.i] < -rhs || t_max_[he.j] < -rhs)) {
      ++n_hold_skip_infeasible;
      continue;
    }
    hold_kept.push_back({he.i, he.j, rhs});
  }

  // Step 2: Fixed-point bound tightening (multi-hop feasibility check)
  // For hold constraint a[j] - a[i] <= rhs:
  //   ub[j] = min(ub[j], ub[i] + rhs)   -- incoming tightens upper bound
  //   lb[i] = max(lb[i], lb[j] - rhs)    -- outgoing tightens lower bound
  // After convergence, prune edges where lb[j] > ub[i] + rhs.
  // This catches multi-hop squeeze infeasibilities that single-edge skip misses.
  // Example: i->j->k where j's ub (from incoming) < j's lb (from outgoing).
  std::vector<float> lb(n_ffs, 0.0f);
  std::vector<float> ub(t_max_.begin(), t_max_.end());

  // Collect negative-rhs edges for propagation (positive rhs cannot tighten)
  std::vector<int> neg_idx;
  neg_idx.reserve(hold_kept.size());
  for (int k = 0; k < static_cast<int>(hold_kept.size()); ++k) {
    if (hold_kept[k].rhs < 0.0f) neg_idx.push_back(k);
  }

  constexpr int BT_MAX_ITER = 50;
  int bt_iter = 0;
  for (; bt_iter < BT_MAX_ITER; ++bt_iter) {
    bool changed = false;
    for (int k : neg_idx) {
      const auto& h = hold_kept[k];
      // ub[j] <= ub[i] + rhs
      float new_ub_j = ub[h.i] + h.rhs;
      if (new_ub_j < ub[h.j] - 1e-6f) {
        ub[h.j] = std::max(0.0f, new_ub_j);
        changed = true;
      }
      // lb[i] >= lb[j] - rhs
      float new_lb_i = lb[h.j] - h.rhs;
      if (new_lb_i > lb[h.i] + 1e-6f) {
        lb[h.i] = std::min(ub[h.i], new_lb_i);
        changed = true;
      }
    }
    if (!changed) break;
  }

  // Clamp squeezed FFs (lb > ub after propagation)
  int n_squeezed = 0;
  for (int idx = 0; idx < n_ffs; ++idx) {
    if (lb[idx] > ub[idx] + 1e-6f) {
      ++n_squeezed;
      float mid = std::max(0.0f, std::min((lb[idx] + ub[idx]) / 2.0f, t_max_[idx]));
      lb[idx] = mid;
      ub[idx] = mid;
    }
  }

  // Prune hold edges using tightened bounds
  int n_bt_pruned = 0;
  for (const auto& h : hold_kept) {
    if (lb[h.j] > ub[h.i] + h.rhs + 1e-6f) {
      ++n_bt_pruned;
      continue;
    }
    MPConstraint* ct = solver->MakeRowConstraint(-inf, static_cast<double>(h.rhs), "");
    ct->SetCoefficient(a_vars[h.j], 1.0);
    ct->SetCoefficient(a_vars[h.i], -1.0);
    ++n_hold_ct;
  }

  if (n_bt_pruned > 0 || n_squeezed > 0) {
    logger_->info(utl::CTS, 658,
        "3D-CTS LP: Bound tightening: {} iters, {} squeezed FFs, {} edges pruned "
        "(from {} neg-rhs edges)",
        bt_iter, n_squeezed, n_bt_pruned, neg_idx.size());
  }

  // Hard PI→FF hold: a_j ≤ corrected_slack  (no penalty variable)
  for (const auto& [ff_idx, slack] : hard_io) {
    MPConstraint* ct = solver->MakeRowConstraint(-inf, slack, "");
    ct->SetCoefficient(a_vars[ff_idx], 1.0);
    ++n_hard_ct;
  }

  // IO soft constraints (with penalty variable w_m)
  for (int m = 0; m < n_io_soft; ++m) {
    const auto& [ff_idx, slack, etype] = soft_io[m];
    MPConstraint* ct = solver->MakeRowConstraint(-inf, slack, "");

    switch (etype) {
      case LpIOEdge::PI_HOLD:
        // a_j - w_m ≤ slack
        ct->SetCoefficient(a_vars[ff_idx], 1.0);
        ct->SetCoefficient(w_vars[m], -1.0);
        break;
      case LpIOEdge::PI_SETUP:
        // -a_j - w_m ≤ slack
        ct->SetCoefficient(a_vars[ff_idx], -1.0);
        ct->SetCoefficient(w_vars[m], -1.0);
        break;
      case LpIOEdge::PO_SETUP:
        // a_i - w_m ≤ slack
        ct->SetCoefficient(a_vars[ff_idx], 1.0);
        ct->SetCoefficient(w_vars[m], -1.0);
        break;
      case LpIOEdge::PO_HOLD:
        // -a_i - w_m ≤ slack
        ct->SetCoefficient(a_vars[ff_idx], -1.0);
        ct->SetCoefficient(w_vars[m], -1.0);
        break;
    }
    ++n_io_ct;
  }

  // Log skeleton stats
  if (n_hold_ghost > 0 || n_hold_skip_infeasible > 0) {
    logger_->info(utl::CTS, CTS_LP_VARS + 2,
        "3D-CTS LP: Skipped {} hold-ghost + {} infeasible hold edges",
        n_hold_ghost, n_hold_skip_infeasible);
  }

  // ===== Solve =====
  auto statusName = [](int s) -> const char* {
    switch (s) {
      case 0: return "OPTIMAL";
      case 1: return "FEASIBLE";
      case 2: return "INFEASIBLE";
      case 3: return "UNBOUNDED";
      case 4: return "ABNORMAL";
      default: return "NOT_SOLVED";
    }
  };

  double elapsed_s = 0.0;
  MPSolver::ResultStatus status = MPSolver::NOT_SOLVED;

  if (enable_cg) {
    // ===== Constraint Generation: setup edges added iteratively =====
    // Hold/IO/hard built upfront. Setup added only when violated.
    // Replaces TopK with exact iterative approach.
    constexpr int CG_MAX_ITER = 30;
    constexpr double CG_VIOL_EPS = 1e-6;

    // Track which pruned_setup edges are already in LP
    std::vector<bool> ps_in_lp(pruned_setup.size(), false);

    // Initial objective: Phase B (TNS minimization)
    // With no setup constraints, a_i = 0 (regularization drives to 0).
    // Hold + IO constraints are the only active ones.
    setPhaseBObjective();

    auto cg_t0 = std::chrono::high_resolution_clock::now();
    int cg_total_added = 0;

    for (int cg_iter = 0; cg_iter < CG_MAX_ITER; ++cg_iter) {
      status = solver->Solve();
      auto cg_t1 = std::chrono::high_resolution_clock::now();
      elapsed_s = std::chrono::duration<double>(cg_t1 - cg_t0).count();

      if (status != MPSolver::OPTIMAL && status != MPSolver::FEASIBLE) {
        logger_->warn(utl::CTS, CTS_LP_RESULT,
            "3D-CTS LP CG iter {}: solver failed (status={})",
            cg_iter, statusName(static_cast<int>(status)));
        break;
      }

      // Cache all arrival values BEFORE modifying the model.
      // Reading solution_value() after model change triggers OR-Tools glog
      // warning "model has been changed" on every call (16M+ lines for ariane133).
      std::vector<double> a_sol(a_vars.size());
      for (size_t i = 0; i < a_vars.size(); ++i) {
        a_sol[i] = a_vars[i]->solution_value();
      }

      // Check ALL pruned_setup edges for violations
      int added = 0;
      for (int idx = 0; idx < static_cast<int>(pruned_setup.size()); ++idx) {
        if (ps_in_lp[idx]) continue;
        const auto& ps = pruned_setup[idx];
        double ai = a_sol[ps.i];
        double aj = a_sol[ps.j];
        if (ai - aj > ps.eff_slack + CG_VIOL_EPS) {
          // V_j: create if new capture endpoint
          if (Vj_map.find(ps.j) == Vj_map.end()) {
            auto* Vj = solver->MakeNumVar(0.0, inf,
                "Vj_" + std::to_string(ps.j));
            Vj_map[ps.j] = Vj;
            // TNS objective
            objective->SetCoefficient(Vj, 1.0);
            // WNS linkage: V_j - W ≤ 0
            if (W_var) {
              MPConstraint* wct = solver->MakeRowConstraint(-inf, 0.0, "");
              wct->SetCoefficient(Vj, 1.0);
              wct->SetCoefficient(W_var, -1.0);
              ++n_wns_ct;
            }
          }

          // Compact setup: a_i - a_j - V_j ≤ eff_slack (no per-edge v_k)
          MPConstraint* ct = solver->MakeRowConstraint(-inf, ps.eff_slack, "");
          ct->SetCoefficient(a_vars[ps.i], 1.0);
          ct->SetCoefficient(a_vars[ps.j], -1.0);
          ct->SetCoefficient(Vj_map.at(ps.j), -1.0);
          ++n_setup_ct;

          ps_in_lp[idx] = true;
          ++added;
        }
      }
      cg_total_added += added;

      const int total_edges = static_cast<int>(pruned_setup.size());
      const double pct = total_edges > 0
          ? 100.0 * cg_total_added / total_edges : 0.0;

      logger_->info(utl::CTS, 702,
          "3D-CTS LP CG iter {}: +{} setup edges (total {}/{} = {:.1f}%, "
          "{} endpoints) in {:.2f}s",
          cg_iter, added, cg_total_added, total_edges, pct,
          static_cast<int>(Vj_map.size()), elapsed_s);

      if (added == 0) break;  // converged — no violated edges

      // Early-exit: if iter 0 added >90% of edges, CG is wasteful — stop
      if (cg_iter == 0 && pct > 90.0) {
        logger_->info(utl::CTS, 706,
            "3D-CTS LP CG: iter 0 added {:.1f}% — early exit (dense graph)", pct);
        break;
      }

      // Tailing prevention: epsilon-violated edges, objective change negligible
      if (added <= 10) {
        logger_->info(utl::CTS, 704,
            "3D-CTS LP CG: tailing ({} edges) �� stopping", added);
        break;
      }
    }

    // After CG convergence: handle 2-phase if enabled
    if (params_.enable_2phase && W_var
        && (status == MPSolver::OPTIMAL || status == MPSolver::FEASIBLE)) {
      // Phase A: minimize W using converged constraint set
      setPhaseAObjective();
      auto ta0 = std::chrono::high_resolution_clock::now();
      MPSolver::ResultStatus st_a = solver->Solve();
      auto ta1 = std::chrono::high_resolution_clock::now();
      elapsed_s += std::chrono::duration<double>(ta1 - ta0).count();

      double W_star = (st_a == MPSolver::OPTIMAL || st_a == MPSolver::FEASIBLE)
                      ? W_var->solution_value() : 0.0;
      logger_->info(utl::CTS, CTS_LP_RESULT,
          "3D-CTS LP CG 2Ph-A: W_star={:.1f}ps in {:.3f}s",
          W_star * 1e3, std::chrono::duration<double>(ta1 - ta0).count());

      // Add WNS bound + Phase B re-solve
      double wns_bound = W_star * (1.0 + static_cast<double>(params_.wns_alpha));
      MPConstraint* wbc = solver->MakeRowConstraint(-inf, wns_bound, "wns_bound");
      wbc->SetCoefficient(W_var, 1.0);

      setPhaseBObjective();
      auto tb0 = std::chrono::high_resolution_clock::now();
      status = solver->Solve();
      auto tb1 = std::chrono::high_resolution_clock::now();
      elapsed_s += std::chrono::duration<double>(tb1 - tb0).count();
    }

  } else {
    // ===== Non-CG: standard solve (legacy path) =====
    int n_ep_vars = static_cast<int>(Vj_map.size());
    int total_ct = n_setup_ct + n_hold_ct + n_hard_ct + n_io_ct + n_wns_ct;
    int total_vars = n_ffs + n_ep_vars + n_io_soft + (W_var ? 1 : 0);

    logger_->info(utl::CTS, CTS_LP_VARS + 1,
        "3D-CTS LP: {} vars (a:{} Vj:{} w:{} W:{}), {} constraints "
        "(setup:{} hold:{} hard_pi:{} io:{} wns:{}) mode={}",
        total_vars, n_ffs, n_ep_vars, n_io_soft, W_var ? 1 : 0,
        total_ct, n_setup_ct, n_hold_ct, n_hard_ct, n_io_ct, n_wns_ct,
        params_.enable_2phase ? "2phase" : "single");

    if (params_.enable_2phase && W_var) {
      setPhaseAObjective();
      auto t0 = std::chrono::high_resolution_clock::now();
      const MPSolver::ResultStatus st_a = solver->Solve();
      auto t1 = std::chrono::high_resolution_clock::now();
      elapsed_s += std::chrono::duration<double>(t1 - t0).count();

      double W_star = (st_a == MPSolver::OPTIMAL || st_a == MPSolver::FEASIBLE)
                      ? W_var->solution_value() : 0.0;
      logger_->info(utl::CTS, CTS_LP_RESULT,
          "3D-CTS LP 2Ph-A: W_star={:.1f}ps in {:.3f}s (status={})",
          W_star * 1e3, std::chrono::duration<double>(t1 - t0).count(),
          static_cast<int>(st_a));

      double wns_bound = W_star * (1.0 + static_cast<double>(params_.wns_alpha));
      MPConstraint* wbc = solver->MakeRowConstraint(-inf, wns_bound, "wns_bound");
      wbc->SetCoefficient(W_var, 1.0);
      logger_->info(utl::CTS, CTS_LP_RESULT,
          "3D-CTS LP 2Ph-B: WNS bound W≤{:.1f}ps (alpha={:.2f})",
          wns_bound * 1e3, params_.wns_alpha);
    }

    setPhaseBObjective();
    auto t_start = std::chrono::high_resolution_clock::now();
    status = solver->Solve();
    auto t_end = std::chrono::high_resolution_clock::now();
    elapsed_s += std::chrono::duration<double>(t_end - t_start).count();
  }

  // ===== Check solver status =====
  result_.solve_time_s = elapsed_s;
  result_.n_ffs = n_ffs;
  result_.n_edges = static_cast<int>(edges_.size());

  if (status != MPSolver::OPTIMAL && status != MPSolver::FEASIBLE) {
    logger_->warn(utl::CTS, CTS_LP_RESULT,
                  "3D-CTS LP: LP solver failed (status={} {})",
                  static_cast<int>(status), statusName(static_cast<int>(status)));
    result_.status = "failed";
    for (const auto& ff : ff_list_) {
      result_.arrivals[ff] = 0.0;
    }
    return result_;
  }

  // ===== Extract results =====
  result_.status = "optimal";
  result_.objective = objective->Value();

  // Per-FF arrivals
  int n_nonzero = 0;
  double sum_arrivals = 0.0;
  double max_arrival = 0.0;
  for (int i = 0; i < n_ffs; ++i) {
    double val = a_vars[i]->solution_value();
    result_.arrivals[ff_list_[i]] = val;
    if (val > 1e-6) {
      ++n_nonzero;
      sum_arrivals += val;
      max_arrival = std::max(max_arrival, val);
    }
  }
  result_.n_nonzero_targets = n_nonzero;

  // Compute endpoint TNS (STA definition) before and after.
  // Use Vj_map (works for both CG and non-CG paths).
  std::unordered_map<int, double> ep_worst_before_map;
  for (const auto& ps : pruned_setup) {
    auto it = ep_worst_before_map.find(ps.j);
    if (it == ep_worst_before_map.end() || ps.eff_slack < it->second) {
      ep_worst_before_map[ps.j] = ps.eff_slack;
    }
  }
  double total_ep_viol_before = 0.0;
  for (auto& [j, s] : ep_worst_before_map) {
    total_ep_viol_before += std::max(0.0, -s);
  }
  double total_ep_viol_after = 0.0;
  for (auto& [j, Vj] : Vj_map) {
    total_ep_viol_after += Vj->solution_value();
  }
  // Compact formulation: no per-edge v_k. total_ep_viol_after already captures TNS.
  double tns_reduction_pct = 0.0;
  if (total_ep_viol_before > 1e-9) {
    tns_reduction_pct = (1.0 - total_ep_viol_after / total_ep_viol_before) * 100.0;
  }

  logger_->info(utl::CTS, CTS_LP_RESULT + 1,
                "3D-CTS LP {}: Solved in {:.2f}s total, obj={:.4f}",
                params_.enable_2phase ? "2Ph-B" : "single",
                elapsed_s, result_.objective);
  logger_->info(utl::CTS, CTS_LP_TARGETS,
                "3D-CTS LP: {} non-zero targets (max={:.1f}ps, sum={:.1f}ps)",
                n_nonzero, max_arrival * 1e3, sum_arrivals * 1e3);
  logger_->info(utl::CTS, CTS_LP_TARGETS,
                "3D-CTS LP: Endpoint TNS (STA): {:.1f}ps → {:.1f}ps ({:.1f}% reduction)",
                total_ep_viol_before * 1e3, total_ep_viol_after * 1e3, tns_reduction_pct);
  // Compact formulation: edge-sum = endpoint TNS (no per-edge v_k)

  // WNS minimax report (endpoint-based)
  if (W_var) {
    double W_val = W_var->solution_value();
    double worst_Vj = 0.0;
    for (auto& [j, Vj] : Vj_map) {
      worst_Vj = std::max(worst_Vj, Vj->solution_value());
    }
    logger_->info(utl::CTS, CTS_LP_TARGETS + 1,
                  "3D-CTS LP: WNS minimax W={:.1f}ps, worst V_j={:.1f}ps (gamma={})",
                  W_val * 1e3, worst_Vj * 1e3, params_.gamma_wns);
  }

  return result_;
}

// Write skew targets CSV for downstream consumers.
// buffer_sizing_lp.py reads pre_cts_skew_targets.csv for skew-aware sizing.
// clock_layer_assignment_v43.tcl reads it for bidirectional layer assignment.
// Without this CSV, both operate with targets=0 → no skew-aware optimization.
// Format matches Python cts_skew_lp.py output: ff_name,target_arrival_ns,tier
void CtsSkewLpSolver::writeTargetsCSV(const std::string& csv_path) const
{
  std::ofstream outfile(csv_path);
  if (!outfile.is_open()) {
    logger_->warn(utl::CTS, CTS_LP_TARGETS + 2,
                  "3D-CTS LP: Cannot write targets CSV: {}", csv_path);
    return;
  }

  outfile << "ff_name,target_arrival_ns,tier\n";
  int n_written = 0;
  int n_nonzero = 0;
  for (int i = 0; i < static_cast<int>(ff_list_.size()); ++i) {
    const std::string& ff_name = ff_list_[i];
    auto it = result_.arrivals.find(ff_name);
    double arrival = (it != result_.arrivals.end()) ? it->second : 0.0;
    int tier = (i < static_cast<int>(ff_tiers_.size())) ? ff_tiers_[i] : 0;
    outfile << ff_name << "," << std::fixed << std::setprecision(6)
            << arrival << "," << tier << "\n";
    ++n_written;
    if (arrival > 1e-6) ++n_nonzero;
  }

  logger_->info(utl::CTS, CTS_LP_TARGETS + 2,
                "3D-CTS LP: Wrote {} targets ({} non-zero) to {}",
                n_written, n_nonzero, csv_path);
}

}  // namespace cts
