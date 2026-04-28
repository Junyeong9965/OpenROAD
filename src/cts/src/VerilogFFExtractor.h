// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors
//
// FF-to-FF Connectivity Extractor
// Two modes: (1) Verilog-based BFS (original), (2) ODB+STA-based (MACRO support)
// BUF_MACRO: Added ODB+STA extraction for macro-included designs.
// Verilog-based BFS cannot see through macro blackboxes; STA-based extraction
// uses findPathEnds which naturally handles macro timing models (Liberty).

#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace odb {
class dbBlock;
class dbInst;
}

namespace sta {
class dbSta;
class dbNetwork;
}

namespace utl {
class Logger;
}

namespace cts {

class Cts3DDatabase;  // Forward decl for 3D tier queries

// IO edge (PI→FF or FF→PO) with worst setup/hold slack.
// Extracted via STA per-port findPathEnds (exhaustive, no worst-N truncation).
struct IOEdgeVerilog {
  std::string edge_type;   // "PI_TO_FF" or "FF_TO_PO"
  std::string port_name;   // Primary input or output port name
  std::string ff_name;     // FF instance name (Verilog format for CSV consistency)
  float slack_setup_ns = 0.0f;  // Worst setup slack (ns)
  float slack_hold_ns = 0.0f;   // Worst hold slack (ns)
  bool has_setup = false;
  bool has_hold = false;
  // Per-clock-net LP: which clock net this FF belongs to.
  std::string ff_clock_net;
};

struct FFEdgeVerilog {
  std::string from_ff;
  std::string to_ff;
  std::string via_net;

  // Locations (will be filled from ODB if available)
  int from_x = 0;
  int from_y = 0;
  int from_tier = 0;  // Tier info for 3D CTS
  int to_x = 0;
  int to_y = 0;
  int to_tier = 0;    // Tier info for 3D CTS

  // Setup timing (max path - will be filled from STA if available)
  float slack_max = 0.0f;
  float arrival_max = 0.0f;
  float required_max = 0.0f;
  // Hold timing (min path)
  float slack_min = 0.0f;
  float arrival_min = 0.0f;
  float required_min = 0.0f;
  bool has_timing = false;  // setup timing filled
  bool has_hold = false;    // hold timing filled
  // Per-clock-net LP: clock net name for each FF endpoint.
  // Populated from ODB (FF clock pin → connected net → trace to root).
  // Used by CtsSkewLpSolver to filter intra-clock vs cross-clock edges.
  std::string from_clock_net;
  std::string to_clock_net;
};

struct RegisterInfo {
  std::string cell_type;
  std::unordered_map<std::string, std::string> pins;  // pin_name -> net_name
};

class VerilogFFExtractor {
public:
  VerilogFFExtractor(const std::string& verilog_path,
                     odb::dbBlock* block,
                     sta::dbSta* sta,
                     sta::dbNetwork* network,
                     utl::Logger* logger);

  // Set 3D database for tier queries (optional)
  void setCts3DDatabase(Cts3DDatabase* db3d) { cts3dDb_ = db3d; }

  // Parse verilog file to extract instances and connections
  void parseVerilog();

  // Extract FF-to-FF edges based on connectivity
  std::vector<FFEdgeVerilog> extractFFEdges();

  // Fill timing info for edges using STA
  void fillTimingInfo(std::vector<FFEdgeVerilog>& edges);
  // Batch timing: capture-side batched findPathEnds (2-5x faster, exact results)
  void fillTimingInfoBatch(std::vector<FFEdgeVerilog>& edges);

  // Write edges to CSV file
  void writeCSV(const std::vector<FFEdgeVerilog>& edges,
                const std::string& output_file);

  // Exhaustive IO edge extraction via STA per-port.
  // PI→FF: for each input port, findPathEnds to all reachable FFs.
  // FF→PO: for each output port, findPathEnds from all reachable FFs.
  // Requires parseVerilog() called first (for register name mapping).
  std::vector<IOEdgeVerilog> extractIOEdges();

  // Write IO edges to CSV (compatible with cts_skew_lp.py parse_io_timing_edges)
  void writeIOCSV(const std::vector<IOEdgeVerilog>& edges,
                  const std::string& output_file);

  // BUF_MACRO: ODB+STA-based FF extraction for macro designs.
  // Collects all sequential instances from ODB (Liberty hasSequentials check),
  // then uses STA findPathEnds to discover all FF→FF timing edges.
  // This replaces parseVerilog()+extractFFEdges()+fillTimingInfo() for macro designs.
  // Advantages over Verilog-based:
  //   - Sees through macro blackboxes (STA uses Liberty timing models)
  //   - No Verilog file needed (works directly from ODB database)
  //   - Handles hierarchical/escaped names correctly (ODB canonical names)
  // The output CSV format is identical to the Verilog-based version.
  void collectRegistersFromODB();
  std::vector<FFEdgeVerilog> extractFFEdgesViaSTA();

private:
  std::string verilog_path_;
  odb::dbBlock* block_;
  sta::dbSta* sta_;
  sta::dbNetwork* network_;
  utl::Logger* logger_;
  Cts3DDatabase* cts3dDb_ = nullptr;  // 3D tier database

  // Parsed data
  std::unordered_map<std::string, RegisterInfo> registers_;
  std::unordered_map<std::string, RegisterInfo> all_instances_;  // All instances including gates
  std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>>
      net_connections_;  // net_name -> [(inst, pin), ...]

  // Helper: Check if cell type is a register
  bool isRegisterCell(const std::string& cell_type) const;

  // Helper: Get location from ODB
  void fillLocations(FFEdgeVerilog& edge);

  // Per-clock-net LP: trace FF's clock pin back to root clock net.
  std::string traceClockNetRoot(odb::dbInst* ff_inst) const;

  // BUF_MACRO: ODB-discovered register instances (name → dbInst*)
  std::unordered_map<std::string, odb::dbInst*> odb_registers_;
};

}  // namespace cts
