// SPDX-License-Identifier: BSD-3-Clause
// C++ LP-based buffer sizing solver
// Replaces Python buffer_sizing_lp.py using OR-Tools GLOP.
//
// Key improvements over Python version:
//   - Liberty-based buffer delay estimation (no analytical 0.024/sqrt(S) model)
//   - Direct STA timing edge access (no CSV round-trip)
//   - Same OR-Tools GLOP solver as CTS LP (proven, fast)
//
// LP formulation:
//   Variables:
//     d_k       = delay of buffer k (continuous, in [d_min_k, d_max_k])
//     p_k, n_k  = positive/negative parts of |d_k - d_k_current| (regularization)
//     s_m, t_m  = positive/negative parts of |d_k - target_k| (skew guidance)
//   Objective (minimize):
//     setup_w * Σ(setup-critical edge weights)
//     + hold_w * Σ(hold-critical edge weights)
//     + reg_w * Σ(p_k + n_k)
//     + skew_w * Σ(s_m + t_m)
//   Constraints:
//     Setup: d_launch - d_capture ≤ setup_slack + margin + current_delta
//     Hold:  d_capture - d_launch ≤ hold_slack - margin + current_delta
//     Linearization: d_k - p_k + n_k ≤ d_current_k (and symmetric)
//     Skew target: d_k - s_m + t_m ≤ target_k (and symmetric)

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace odb {
class dbBlock;
class dbMaster;
}

namespace sta {
class dbSta;
class dbNetwork;
class LibertyCell;
class LibertyPort;
class TimingArc;
class Corner;
}

namespace utl {
class Logger;
}

namespace cts {

// Buffer info extracted from ODB
struct BufInfo {
  std::string inst_name;
  std::string cell_name;
  int buf_idx;
  std::vector<std::string> driven_ffs;  // FF names driven by this buffer
  double current_delay_ps;              // Liberty-measured delay
};

// Sizing candidate: (cell_name, delay_ps)
struct SizeOption {
  std::string cell_name;
  double delay_ps;
};

// Sizing change result
struct SizingChange {
  std::string inst_name;
  std::string old_cell;
  std::string new_cell;
  std::string direction;  // "upsize" or "downsize"
};

// Buffer sizing LP parameters
struct BufSizingParams {
  double setup_margin_ps = 0.0;
  double hold_margin_ps = 0.0;
  double hold_weight = 0.5;
  double reg_weight = 0.01;
  double skew_weight = 0.5;
};

class BufSizingLpSolver {
 public:
  BufSizingLpSolver(odb::dbBlock* block,
                    sta::dbSta* sta,
                    utl::Logger* logger);

  void setParams(const BufSizingParams& params) { params_ = params; }

  // Step 1: Extract clock buffers from ODB (matching *clkbuf* pattern)
  void extractBuffers();

  // Step 2: Extract timing edges from STA (setup + hold slacks)
  void extractTimingEdges(int max_paths = 5000);

  // Step 3: Measure per-master Liberty delay (buf driving own input cap)
  void measureLibraryDelays();

  // Step 4: Load LP skew targets (optional, from post-CTS LP CSV)
  void loadSkewTargets(const std::string& csv_path);

  // Step 5: Solve LP
  // Returns sizing changes sorted by priority (largest slack improvement first)
  std::vector<SizingChange> solve();

  // Write sizing results to CSV (for Tcl iterative application)
  void writeResultsCSV(const std::string& csv_path) const;

  // Get measured delay for a buffer master (for external use)
  double getMasterDelay(const std::string& master_name) const;

 private:
  // Measure single Liberty cell delay at given load cap
  double measureCellDelay(sta::LibertyCell* cell, double load_cap);

  // Get all size options for a buffer by tier suffix (_upper/_bottom).
  // Uses tier-based matching (not regex family matching) to support
  // cross-family sizing in heterogeneous 3D (e.g., BUF_X4_upper -> BUFx2_ASAP7...).
  std::vector<SizeOption> getSizeOptions(const std::string& cell_name);

  // Get drive strength from cell name (e.g., "BUFx8_ASAP7..." -> 8)
  static int getDriveStrength(const std::string& name);

  odb::dbBlock* block_;
  sta::dbSta* sta_;
  utl::Logger* logger_;

  BufSizingParams params_;

  // Buffer data
  std::vector<BufInfo> buffers_;
  std::unordered_map<std::string, int> ff_to_buf_;  // ff_name -> buf_idx

  // Timing edges
  struct TimingEdge {
    std::string launch_ff;
    std::string capture_ff;
    double slack_setup_ps;
    double slack_hold_ps;
  };
  std::vector<TimingEdge> edges_;

  // Liberty-measured delays: master_name -> delay_ps
  std::unordered_map<std::string, double> master_delays_;

  // Per-buffer skew targets (from post-CTS LP)
  std::unordered_map<int, double> buf_skew_targets_;  // buf_idx -> target_ps

  // Sizing results
  std::vector<SizingChange> results_;
};

}  // namespace cts
