// SPDX-License-Identifier: BSD-3-Clause
// C++ LP solver for 3D CTS skew targeting
// Endpoint-TNS objective (STA TNS aligned)
// Replaces Python cts_skew_lp.py using Google OR-Tools GLOP
//
// LP-TNS formulation (endpoint-TNS, aligned with STA TNS definition):
//   Variables: a_i (per-FF arrival), v_k (per-edge violation, tie-break only),
//              V_j (per-endpoint worst violation = STA TNS contribution),
//              w_m (IO penalty), W (WNS minimax)
//   Single-phase objective:
//     min lambda_reg*Σa_i + Σ V_j + ε·Σv_k + weight_io*Σw_m + gamma_wns*W
//     Σ V_j aligns LP with STA TNS = Σ_j max(0, WNS_j).
//     ε = 0.001: tie-breaking so that among solutions with same Σ V_j, smaller v_k preferred.
//   2-phase objective (enable_2phase=true):
//     Phase A: min W  (find best achievable WNS → W_star)
//     Phase B: min Σ V_j + ε·Σv_k + lambda_reg*Σa_i + weight_io*Σw_m
//              s.t. W ≤ W_star*(1+wns_alpha)  (WNS bound, default alpha=0.05)
//     Advantage: no gamma_wns hyperparameter needed; explicit WNS guarantee.
//   Constraints (both modes):
//     FF→FF setup: a_i - a_j - v_k ≤ eff_setup_k     (soft via v_k)
//     FF→FF hold:  a_j - a_i ≤ eff_hold_k - margin    (HARD)
//     PI→FF hold:  a_j ≤ corrected_slack              (HARD, positive slack only)
//     IO soft:     ±a_j - w_m ≤ slack                  (soft via w_m)
//     Endpoint:    v_k - V_j ≤ 0  for all k arriving at j   (forces V_j = max v_k→j)
//     WNS minimax: V_j - W ≤ 0                        (endpoint-based, N_ep constraints)
//   Bounds: 0 ≤ a_i ≤ t_max_i

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace odb {
class dbBlock;
}

namespace sta {
class dbSta;
class dbNetwork;
}

namespace utl {
class Logger;
}

namespace cts {

class Cts3DDatabase;
class VerilogFFExtractor;

// Timing edge between two FFs (from STA)
struct LpTimingEdge {
  int from_idx;
  int to_idx;
  float slack_setup_ns;  // max path slack (setup)
  float slack_hold_ns;   // min path slack (hold)
  bool cross_tier;
  bool hold_ghost;       // true if STA found no hold path (slack_min=0, arrival_min=0)
};

// IO edge (PI→FF or FF→PO)
struct LpIOEdge {
  int ff_idx;
  float slack_ns;
  enum Type { PI_HOLD, PI_SETUP, PO_SETUP, PO_HOLD } type;
};

// LP solver parameters
struct LpParams {
  float sigma_local = 0.005f;    // local clock uncertainty (ns)
  float sigma_pi = 0.005f;      // PI→FF hold safety margin (ns)
  float lambda_reg = 0.01f;     // regularization weight (matched to Python default)
  float hold_margin = 0.010f;   // hold margin (ns)
  float gamma_wns = 10.0f;      // WNS minimax penalty weight (Phase B tie-break)
  float weight_io = 1.0f;       // IO soft constraint weight
  float max_skew = 0.100f;      // max achievable skew (ns)
  bool hard_pi_hold = false;    // use hard PI→FF hold constraints
  // 2-phase LP
  // Phase A: minimize W (best achievable WNS) → W_star
  // Phase B: minimize Σ V_j + ε·Σv_k s.t. W ≤ W_star*(1+wns_alpha)
  bool enable_2phase = false;   // enable 2-phase LP (false = single-phase)
  float wns_alpha = 0.05f;      // WNS relaxation factor in Phase B bound
  // Delivery bound tightening
  // Clamp t_max to MAX_DEPTH * d_buf — prevents LP from assigning
  // targets beyond cascaded TAP physical delivery capacity.
  // 0 = disabled (use max_skew only).
  float delivery_max = 0.0f;
  // Endpoint-TopK pruning: per capture endpoint j, keep only the worst K
  // incoming setup edges (by eff_slack).  Directly controls LP setup constraint
  // count: n_capture_endpoints * K.  Hold edges are NOT pruned.
  // 0 = disabled (keep all setup edges).  Default 16.
  int endpoint_topk = 0;  // 0=disabled (keep all). Set via CTS_LP_ENDPOINT_TOPK for large designs.
};

// LP solver result
struct LpResult {
  std::unordered_map<std::string, double> arrivals;  // ff_name -> arrival_ns
  std::string status;
  double objective = 0.0;
  double solve_time_s = 0.0;
  int n_ffs = 0;
  int n_edges = 0;
  int n_nonzero_targets = 0;
};

class CtsSkewLpSolver {
 public:
  CtsSkewLpSolver(odb::dbBlock* block,
                  sta::dbSta* sta,
                  sta::dbNetwork* network,
                  utl::Logger* logger,
                  Cts3DDatabase* cts3dDb);

  // Set LP parameters
  void setParams(const LpParams& params) { params_ = params; }

  // Read FF-to-FF timing graph from Phase 1 CSV.
  // clock_net_filter: if non-empty, only include edges where BOTH from and to
  // belong to this clock net. Cross-clock edges are excluded.
  void readTimingGraphCSV(const std::string& csv_path,
                          const std::string& clock_net_filter = "");

  // Phase 2: Read IO edges from Phase 1 CSV (PI→FF, FF→PO)
  void readIOEdgesCSV(const std::string& io_csv_path,
                      const std::string& clock_net_filter = "");

  // Phase 3: Compute per-FF bounds from clock latency estimator
  void computeBounds();

  // Phase 4: Apply PI-hold-clip (reduce t_max for PI→FF hold-critical FFs)
  void applyPiHoldClip();

  // Phase 5: Solve LP-TNS
  LpResult solveTns();

  // Get results for direct consumption by HTreeBuilder
  const std::unordered_map<std::string, double>& getArrivals() const {
    return result_.arrivals;
  }

  // Get number of FFs
  int getNumFFs() const { return static_cast<int>(ff_list_.size()); }

  // Write skew targets to CSV (for downstream buffer_sizing_lp.py and layer assignment)
  void writeTargetsCSV(const std::string& csv_path) const;

 private:
  odb::dbBlock* block_;
  sta::dbSta* sta_;
  sta::dbNetwork* network_;
  utl::Logger* logger_;
  Cts3DDatabase* cts3dDb_;

  LpParams params_;
  LpResult result_;

  // Verilog path for VerilogFFExtractor reuse
  std::string verilog_path_;

  // FF data
  std::vector<std::string> ff_list_;       // ordered FF names
  std::unordered_map<std::string, int> ff_to_idx_;
  std::vector<int> ff_tiers_;              // tier per FF
  std::vector<float> t_max_;              // per-FF max achievable arrival (ns)

  // Timing edges (in-memory, no CSV)
  std::vector<LpTimingEdge> edges_;
  std::vector<LpIOEdge> io_edges_;

  // Ghost edge filter: skip edges with all-zero timing
  bool isGhostEdge(float slack_max, float slack_min,
                   float arr_max, float arr_min) const;
};

}  // namespace cts
