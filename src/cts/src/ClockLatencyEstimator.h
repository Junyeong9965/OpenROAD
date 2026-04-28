// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors
//
// ClockLatencyEstimator — per-FF physical achievability
// bounds estimator for the LP-SAFETY pre-CTS skew solver.
//
// Implements the professor's "three leaf buffer cases" (bottom-*):
//   Case 1: Single bottom-tier leaf buffer drives:
//             - bottom FFs directly (delay in [0, max_skew])
//             - top    FFs through HB via (delay in [t_via, t_via + max_skew])
//   Case 2: Single top-tier leaf buffer drives:
//             - top    FFs directly (delay in [0, max_skew])
//             - bottom FFs through HB via (delay in [t_via, t_via + max_skew])
//   Case 3: Per-tier leaf buffers (top + bottom independent):
//             - each FF driven on its own tier (delay in [0, max_skew])
//
// Union of achievable delay ranges across all three cases:
//   For any FF on either tier: a_i in [0, max_skew + t_via]
//   where t_via = 0.693 * R_HB * C_HB (Elmore 50% delay of HB via)
//
// Output CSV format: ff_name,tier,t_min_ns,t_max_ns
//   t_min_ns = 0.0        for all FFs  (add-only delay; a_i >= 0)
//   t_max_ns = max_skew   for bottom-tier FFs
//   t_max_ns = max_skew + t_via  for upper-tier FFs

#pragma once

#include <string>

namespace cts {

class Cts3DDatabase;

class ClockLatencyEstimator
{
 public:
  // Construct with a populated Cts3DDatabase (tiers + HBT parasitic loaded).
  explicit ClockLatencyEstimator(Cts3DDatabase* db3d);

  // Compute per-FF physical achievability bounds and write CSV.
  // max_skew_ns: CTS delay budget per FF for same-tier tree sizing (ns).
  // output_csv:  destination path, format ff_name,tier,t_min_ns,t_max_ns.
  void estimateAndWrite(const std::string& output_csv, double max_skew_ns);

  // HB via delay in ns (0.693 * R_HB * C_HB, computed at construction).
  double getViaDelayNs() const { return t_via_ns_; }

 private:
  Cts3DDatabase* db3d_;    // 3D database: tier assignments + HBT R/C
  double t_via_ns_ = 0.0;  // Elmore delay of one HB via (ns)
};

}  // namespace cts
