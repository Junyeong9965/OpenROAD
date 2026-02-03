// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors
//
// Verilog-based FF-to-FF Connectivity Extractor
// Parse verilog netlist to extract FF-to-FF edges with timing info from STA

#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
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

struct FFEdgeVerilog {
  std::string from_ff;
  std::string to_ff;
  std::string via_net;

  // Locations (will be filled from ODB if available)
  int from_x = 0;
  int from_y = 0;
  int from_tier = 0;
  int to_x = 0;
  int to_y = 0;
  int to_tier = 0;

  // Setup timing (max path - will be filled from STA if available)
  float slack_max = 0.0f;
  float arrival_max = 0.0f;
  float required_max = 0.0f;
  // Hold timing (min path)
  float slack_min = 0.0f;
  float arrival_min = 0.0f;
  float required_min = 0.0f;
  bool has_timing = false;
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

  // Parse verilog file to extract instances and connections
  void parseVerilog();

  // Extract FF-to-FF edges based on connectivity
  std::vector<FFEdgeVerilog> extractFFEdges();

  // Fill timing info for edges using STA
  void fillTimingInfo(std::vector<FFEdgeVerilog>& edges);

  // Write edges to CSV file
  void writeCSV(const std::vector<FFEdgeVerilog>& edges,
                const std::string& output_file);

private:
  std::string verilog_path_;
  odb::dbBlock* block_;
  sta::dbSta* sta_;
  sta::dbNetwork* network_;
  utl::Logger* logger_;

  // Parsed data
  std::unordered_map<std::string, RegisterInfo> registers_;
  std::unordered_map<std::string, RegisterInfo> all_instances_;  // All instances including gates
  std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>>
      net_connections_;  // net_name -> [(inst, pin), ...]

  // Helper: Check if cell type is a register
  bool isRegisterCell(const std::string& cell_type) const;

  // Helper: Get location from ODB
  void fillLocations(FFEdgeVerilog& edge);
};

}  // namespace cts
