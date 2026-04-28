// SPDX-License-Identifier: BSD-3-Clause
// C++ buffer sizing LP using OR-Tools GLOP
// Replaces Python buffer_sizing_lp.py
//
// Post-CTS buffer sizing: adjusts clock buffer delays to improve timing.
// Variables: d_k (delay per buffer, continuous in [d_min, d_max])
// Objective: min TNS (weighted setup+hold improvement)
// Constraints: setup/hold slack preservation per timing edge
// Post-LP: map continuous delay to nearest discrete buffer size

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace odb {
class dbBlock;
class dbInst;
class dbMaster;
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

// Buffer size option with delay and capacitance
struct BufSizeOption {
  std::string master_name;
  float delay_ps;         // intrinsic delay (ps)
  int drive_strength;
  float input_cap_ff;     // estimated input cap (fF)
};

// Clock buffer info (extracted from ODB after CTS)
struct ClockBufInfo {
  std::string inst_name;
  std::string cell_name;
  std::string tier;  // "bottom" or "upper"
  std::vector<std::string> driven_ffs;
  float current_delay_ps;
  std::vector<BufSizeOption> size_options;  // Pareto-optimal sizes
};

// Buffer sizing result: old → new master
struct BufSizingChange {
  std::string inst_name;
  std::string old_master;
  std::string new_master;
  float old_delay_ps;
  float new_delay_ps;
};

class CtsBufferSizingLp {
 public:
  CtsBufferSizingLp(odb::dbBlock* block,
                    sta::dbSta* sta,
                    sta::dbNetwork* network,
                    utl::Logger* logger);

  // Extract post-CTS timing edges from STA
  void extractTimingEdges();

  // Extract clock buffer inventory from ODB
  void extractBufferInfo();

  // Enumerate available buffer sizes from Liberty
  void enumerateBufferSizes();

  // Solve LP and return sizing changes
  std::vector<BufSizingChange> solve(float hold_weight = 0.5f,
                                     float reg_weight = 0.01f);

  // Apply changes directly to ODB (swapMaster)
  void applyChanges(const std::vector<BufSizingChange>& changes);

 private:
  odb::dbBlock* block_;
  sta::dbSta* sta_;
  sta::dbNetwork* network_;
  utl::Logger* logger_;

  std::vector<ClockBufInfo> buffers_;
  std::unordered_map<std::string, int> ff_to_buf_;  // FF → buffer index

  struct TimingEdge {
    std::string from_ff;
    std::string to_ff;
    float slack_setup_ps;
    float slack_hold_ps;
  };
  std::vector<TimingEdge> edges_;

  // Helper: get buffer delay from Liberty
  float getBufferDelay(const std::string& master_name) const;

  // Helper: is this instance a clock buffer?
  bool isClockBuffer(odb::dbInst* inst) const;
};

}  // namespace cts
