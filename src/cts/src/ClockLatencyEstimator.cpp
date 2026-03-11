// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors
//
// JYJ (2026-02-23) V32: ClockLatencyEstimator implementation.
// See ClockLatencyEstimator.h for design rationale.

#include "ClockLatencyEstimator.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>

#include "Cts3DDatabase.h"
#include "odb/db.h"

namespace cts {

ClockLatencyEstimator::ClockLatencyEstimator(Cts3DDatabase* db3d)
    : db3d_(db3d)
{
  // Compute HB via delay using Elmore delay model (lumped RC, 50% threshold):
  //   t_via = 0.693 * R_HB (ohms) * C_HB (farads)  [seconds]
  // R_HB and C_HB come from set_layer_rc -via hb_layer in the tech setup.
  // If HBT parasitics are not set, t_via = 0 and upper-tier FFs get the
  // same bound as bottom-tier FFs (conservative fallback).
  const double R = db3d_->getHbtResistance();  // ohms
  const double C = db3d_->getHbtCapacitance();  // farads
  t_via_ns_ = 0.693 * R * C * 1.0e9;           // convert seconds -> ns
}

void ClockLatencyEstimator::estimateAndWrite(const std::string& output_csv,
                                              double max_skew_ns)
{
  // Compute t_max for each tier.
  // Union of all three leaf buffer cases (see header for details):
  //   bottom-tier FF: t_max = max_skew   (Case 3 or direct path in Case 1/2)
  //   upper-tier  FF: t_max = max_skew + t_via  (via HB path in Case 1 or 2)
  // t_min = 0 for all FFs (add-only delay model: a_i >= 0 in LP)
  const double t_max_bottom = max_skew_ns;
  const double t_max_upper  = max_skew_ns + t_via_ns_;

  std::ofstream f(output_csv);
  if (!f.is_open()) {
    // Caller (TritonCTS::estimateLeafLatencies) will warn about missing file
    return;
  }

  f << std::fixed << std::setprecision(6);
  f << "ff_name,tier,t_min_ns,t_max_ns\n";

  // Bottom-tier instances (tier 0)
  const auto& tier0 = db3d_->getInstancesOnTier(0);
  for (auto* inst : tier0) {
    f << inst->getName() << ",0,0.000000," << t_max_bottom << "\n";
  }

  // Upper-tier instances (tier 1)
  const auto& tier1 = db3d_->getInstancesOnTier(1);
  for (auto* inst : tier1) {
    f << inst->getName() << ",1,0.000000," << t_max_upper << "\n";
  }
}

}  // namespace cts
