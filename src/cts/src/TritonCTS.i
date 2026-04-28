// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

%{
#include <cstdint>

#include "cts/TritonCTS.h"
#include "CtsOptions.h"
#include "TechChar.h"
#include "CtsGraphics.h"
#include "ord/OpenRoad.hh"

namespace ord {
// Defined in OpenRoad.i
cts::TritonCTS *
getTritonCts();
}

using namespace cts;
using ord::getTritonCts;
%}

%include "../../Exception.i"
%include "stdint.i"
%include "std_string.i"

%ignore cts::CtsOptions::setObserver;
%ignore cts::CtsOptions::getObserver;

// Enum: CtsOptions::NdrStrategy
%typemap(typecheck) CtsOptions::NdrStrategy {
  char *str = Tcl_GetStringFromObj($input, 0);
  if (strcasecmp(str, "NONE") == 0) {
    $1 = 1;
  } else if (strcasecmp(str, "ROOT_ONLY") == 0) {
    $1 = 1;
  } else if (strcasecmp(str, "HALF") == 0) {
    $1 = 1;
  } else if (strcasecmp(str, "FULL") == 0) {
    $1 = 1;
  } else {
    $1 = 0;
  }
}

%typemap(in) CtsOptions::NdrStrategy {
  char *str = Tcl_GetStringFromObj($input, 0);
  if (strcasecmp(str, "ROOT_ONLY") == 0) {
    $1 = CtsOptions::NdrStrategy::ROOT_ONLY;
  } else if (strcasecmp(str, "HALF") == 0) {
    $1 = CtsOptions::NdrStrategy::HALF;
  } else if (strcasecmp(str, "FULL") == 0) {
    $1 = CtsOptions::NdrStrategy::FULL;
  } else {
    $1 = CtsOptions::NdrStrategy::NONE;
  };
}

%inline %{

void
set_sink_clustering(bool enable)
{
  getTritonCts()->getParms()->setSinkClustering(enable);
}

void
set_plot_option(bool plot)
{
  getTritonCts()->getParms()->setPlotSolution(plot);
}

void
set_sink_clustering_levels(unsigned levels)
{
  getTritonCts()->getParms()->setSinkClusteringLevels(levels);
}

void
set_max_char_cap(double cap)
{
  getTritonCts()->getParms()->setMaxCharCap(cap);
}

void
set_max_char_slew(double slew)
{
  getTritonCts()->getParms()->setMaxCharSlew(slew);
}

void
set_wire_segment_distance_unit(unsigned unit)
{
  getTritonCts()->getParms()->setWireSegmentUnit(unit);
}

void
set_root_buffer(const char* buffer)
{
  getTritonCts()->setRootBuffer(buffer);
}

void
set_slew_steps(int steps)
{
  getTritonCts()->getParms()->setSlewSteps(steps);
}

void
set_cap_steps(int steps)
{
  getTritonCts()->getParms()->setCapSteps(steps);
}

void
set_metric_output(const char* file)
{
  getTritonCts()->getParms()->setMetricsFile(file);
}

void
set_debug_cmd()
{
  getTritonCts()->getParms()->setObserver(std::make_unique<CtsGraphics>());
}

void
report_cts_metrics()
{
  getTritonCts()->reportCtsMetrics();
}

void
set_tree_buf(const char* buffer)
{
  getTritonCts()->getParms()->setTreeBuffer(buffer);
}

void
set_distance_between_buffers(int distance)
{
  getTritonCts()->getParms()->setSimpleSegmentsEnabled(true);
  getTritonCts()->getParms()->setBufferDistance(distance);
}

void
set_branching_point_buffers_distance(int distance)
{
  getTritonCts()->getParms()->setVertexBuffersEnabled(true);
  getTritonCts()->getParms()->setVertexBufferDistance(distance);
}

void
set_clustering_exponent(unsigned power)
{
  getTritonCts()->getParms()->setClusteringPower(power);
}

void
set_clustering_unbalance_ratio(double ratio)
{
  getTritonCts()->getParms()->setClusteringCapacity(ratio);
}

void
set_sink_clustering_size(unsigned size)
{
  getTritonCts()->getParms()->setSinkClusteringSize(size);
}

void
set_clustering_diameter(double distance)
{
  getTritonCts()->getParms()->setMaxDiameter(distance);
}

void
set_macro_clustering_size(unsigned size)
{
  getTritonCts()->getParms()->setMacroClusteringSize(size);
}

void
set_macro_clustering_diameter(double distance)
{
  getTritonCts()->getParms()->setMacroMaxDiameter(distance);
}

void
set_num_static_layers(unsigned num)
{
  getTritonCts()->getParms()->setNumStaticLayers(num);
}

void
set_sink_buffer(const char* buffer)
{
  getTritonCts()->setSinkBuffer(buffer);
}

void
report_characterization()
{
  getTritonCts()->getCharacterization()->report();
}

void
report_wire_segments(unsigned length, unsigned load, unsigned outputSlew)
{
  getTritonCts()->getCharacterization()->reportSegments(length, load, outputSlew);
}

int
set_clock_nets(const char* names)
{
  return getTritonCts()->setClockNets(names);
}

void
set_skip_clock_nets(odb::dbNet* net)
{
  getTritonCts()->getParms()->setSkipNets(net);
}

void
set_buffer_list(const char* buffers)
{
  getTritonCts()->setBufferList(buffers);
}

void
set_obstruction_aware(bool obs)
{
  getTritonCts()->getParms()->setObstructionAware(obs);
}

void
set_apply_ndr(CtsOptions::NdrStrategy strategy)
{
  getTritonCts()->getParms()->setApplyNDR(strategy);
}

void
set_insertion_delay(bool insDelay)
{
  getTritonCts()->getParms()->enableInsertionDelay(insDelay);
}

void
set_sink_buffer_max_cap_derate(double derate)
{
  getTritonCts()->getParms()->setSinkBufferMaxCapDerate(derate);
}

void
set_dummy_load(bool dummyLoad)
{
  getTritonCts()->getParms()->enableDummyLoad(dummyLoad);
}

void
set_delay_buffer_derate(float derate)
{
  getTritonCts()->getParms()->setDelayBufferDerate(derate);
}

void
set_cts_library(const char* name)
{
  getTritonCts()->getParms()->setCtsLibrary(name);
}

void
set_repair_clock_nets(bool value)
{
  getTritonCts()->getParms()->setRepairClockNets(value);
}
 
void
run_triton_cts()
{
  getTritonCts()->runTritonCts();
}

const char*
get_ndr_strategy()
{
  return getTritonCts()->getParms()->getApplyNdrName();
}

int
get_branching_buffers_distance()
{
  return getTritonCts()->getParms()->getVertexBufferDistance();
}
std::string
get_buffer_list()
{
  return getTritonCts()->getParms()->getBufferListToString();
}
unsigned
get_clustering_exponent()
{
  return getTritonCts()->getParms()->getClusteringPower();
}

double
get_clustering_unbalance_ratio()
{
  return getTritonCts()->getParms()->getClusteringCapacity();
}

float
get_delay_buffer_derate()
{
  return getTritonCts()->getParms()->getDelayBufferDerate();
}

int
get_distance_between_buffers()
{
  return getTritonCts()->getParms()->getBufferDistance();
}

const char*
get_library()
{
  return getTritonCts()->getParms()->getCtsLibrary();
}

double
get_macro_clustering_max_diameter()
{
  return getTritonCts()->getParms()->getMacroMaxDiameter();
}

unsigned
get_macro_clustering_size()
{
  return getTritonCts()->getParms()->getMacroSinkClusteringSize();
}

unsigned
get_num_static_layers()
{
  return getTritonCts()->getParms()->getNumStaticLayers();
}

std::string
get_root_buffer()
{
  return getTritonCts()->getRootBufferToString();
}

double
get_sink_buffer_max_cap_derate()
{
  return getTritonCts()->getParms()->getSinkBufferMaxCapDerate();
}

unsigned
get_sink_clustering_levels()
{
  return getTritonCts()->getParms()->getSinkClusteringLevels();
}

double
get_sink_clustering_max_diameter()
{
  return getTritonCts()->getParms()->getMaxDiameter();
}

unsigned
get_sink_clustering_size()
{
  return getTritonCts()->getParms()->getSinkClusteringSize();
}

std::string
get_skip_nets()
{
  return getTritonCts()->getParms()->getSkipNetsToString();
}
std::string
get_tree_buf()
{
  return getTritonCts()->getParms()->getTreeBuffer();
}

unsigned
get_wire_unit()
{
  return getTritonCts()->getParms()->getWireSegmentUnit();
}

void
set_db_unit(bool reset)
{
  if (reset) {
    getTritonCts()->getParms()->setDbUnits(-1);
  } else {
    odb::dbBlock* block = getTritonCts()->getBlock();
    getTritonCts()->getParms()->setDbUnits(block->getDbUnitsPerMicron());
  }
}
void
reset_apply_ndr()
{
  getTritonCts()->getParms()->resetApplyNDR();
}
void
reset_buffer_list()
{
  getTritonCts()->getParms()->resetBufferList();
}
void
reset_branching_point_buffers_distance()
{
  getTritonCts()->getParms()->resetVertexBufferDistance();
}
void
reset_clustering_exponent()
{
  getTritonCts()->getParms()->resetClusteringPower();
}
void
reset_clustering_unbalance_ratio()
{
  getTritonCts()->getParms()->resetClusteringCapacity();
}
void
reset_delay_buffer_derate()
{
  getTritonCts()->getParms()->resetDelayBufferDerate();
}
void
reset_distance_between_buffers()
{
  getTritonCts()->getParms()->resetBufferDistance();
}
void
reset_cts_library()
{
  getTritonCts()->getParms()->resetCtsLibrary();
}
void
reset_macro_clustering_diameter()
{
  getTritonCts()->getParms()->resetMacroMaxDiameter();
}
void
reset_macro_clustering_size()
{
  getTritonCts()->getParms()->resetMacroClusteringSize();
}
void
reset_num_static_layers()
{
  getTritonCts()->getParms()->resetNumStaticLayers();
}
void
reset_root_buffer()
{
  getTritonCts()->resetRootBuffer();
}
void
reset_sink_buffer_max_cap_derate()
{
  getTritonCts()->getParms()->resetSinkBufferMaxCapDerate();
}
void
reset_sink_clustering_levels()
{
  getTritonCts()->getParms()->resetSinkClusteringLevels();
}
void
reset_clustering_diameter()
{
  getTritonCts()->getParms()->resetMaxDiameter();
}
void
reset_sink_clustering_size()
{
  getTritonCts()->getParms()->resetSinkClusteringSize();
}
void
reset_skip_nets()
{
  getTritonCts()->getParms()->resetSkipNets();
}
void
reset_tree_buf()
{
  getTritonCts()->getParms()->resetTreeBuffer();
}
void
reset_wire_segment_distance_unit()
{
  getTritonCts()->getParms()->resetWireSegmentUnit();
}

void
extract_ff_timing_graph_verilog(const char* verilog_file, const char* output_file)
{
  getTritonCts()->extractFFGraphFromVerilog(verilog_file, output_file);
}

// ODB+STA-based FF timing graph extraction.
// No Verilog file needed — registers discovered from ODB Liberty cells,
// edges extracted via STA findPathEnds. Handles macro-included designs.
void
extract_ff_timing_graph_odb(const char* output_file)
{
  getTritonCts()->extractFFGraphFromODB(output_file);
}

// Exhaustive IO timing edge extraction via STA per-port.
// Replaces Tcl extract_io_timing_edges_phase1a (worst-N truncation).
void
extract_io_timing_edges(const char* verilog_file, const char* output_file)
{
  getTritonCts()->extractIOTimingEdges(verilog_file, output_file);
}

// ODB-based IO timing edge extraction.
// No Verilog file needed — uses ODB register collection + STA per-port.
// Handles macro-included designs where Verilog parsing fails.
void
extract_io_timing_edges_odb(const char* output_file)
{
  getTritonCts()->extractIOTimingEdgesFromODB(output_file);
}

// Pre-CTS skew targets for SG-CTS
void
load_skew_targets(const char* csv_file)
{
  getTritonCts()->loadSkewTargets(csv_file);
}

// Estimate per-FF physical achievability bounds.
// Computes HB via delay (t_via = 0.693 * R_HB * C_HB) from tech and writes
// bounds CSV for the LP-SAFETY pre-CTS skew solver (pre_cts_skew_lp.py).
// output_csv:  destination path (ff_name,tier,t_min_ns,t_max_ns)
// max_skew_ns: same-tier CTS delay budget in ns (default 0.1 = 100ps)
void
estimate_leaf_latencies(const char* output_csv, float max_skew_ns)
{
  getTritonCts()->estimateLeafLatencies(output_csv, max_skew_ns);
}

// C++ LP solver for skew targeting.
// Extracts timing graph from STA in-memory, solves LP-TNS with OR-Tools GLOP,
// stores results in Cts3DDatabase skewTargetMap_ (no CSV round-trip).
void
solve_skew_lp(const char* verilog_file,
              float sigma_local, float sigma_pi,
              float lambda_reg, float hold_margin,
              float gamma_wns, float weight_io,
              int hard_pi_hold, float max_skew)
{
  getTritonCts()->solveSkewLp(verilog_file,
                              sigma_local, sigma_pi,
                              lambda_reg, hold_margin,
                              gamma_wns, weight_io,
                              hard_pi_hold != 0, max_skew);
}

// Old solve_buffer_sizing_lp(float, float) removed — replaced by V53_FM_BUF version below.

// Update cascaded TAP chain connections in ODB.
// Re-reads LP targets CSV and reconnects FFs to the correct TAP depth level
// without rebuilding the tree.
void
update_tap_depths(const char* targets_csv)
{
  getTritonCts()->updateTapDepths(targets_csv);
}

// C++ LP-based buffer sizing with Liberty delays.
// Replaces Python buffer_sizing_lp.py. Outputs sizing decisions CSV.
void
solve_buffer_sizing_lp(const char* timing_csv,
                       const char* output_csv,
                       const char* skew_targets_csv,
                       double hold_weight,
                       double reg_weight,
                       double skew_weight,
                       double setup_margin_ps,
                       double hold_margin_ps)
{
  // Set timing CSV path as env var for BufSizingLpSolver to read
  setenv("BUF_SIZING_TIMING_CSV", timing_csv, 1);
  getTritonCts()->solveBufferSizingLp(output_csv, skew_targets_csv,
                                       hold_weight, reg_weight, skew_weight,
                                       setup_margin_ps, hold_margin_ps);
}

%} //inline
