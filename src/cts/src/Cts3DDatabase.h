// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors
//
// JYJ (2026-02-06) Created Cts3DDatabase: centralized 3D tier management
// for True 3D CTS. Replaces scattered tier logic in HTreeBuilder and
// VerilogFFExtractor with a single source of truth.
// JYJ (2026-02-23) V32a: Added explicit tier buffer pair mapping to fix
// getBufferForTier() fallback when cell names differ across technology nodes
// (e.g., "BUF_X4_bottom" -> "BUFx4_ASAP7_75t_R_upper" instead of
// the non-existent "BUF_X4_upper").

#pragma once

#include <array>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "Clock.h"
#include "Util.h"
#include "odb/db.h"
#include "utl/Logger.h"

namespace cts {

class Cts3DDatabase
{
 public:
  Cts3DDatabase(odb::dbDatabase* db, utl::Logger* logger);

  // Load tier info from ODB (call once after design is loaded)
  // Scans all instances, determines tier from master name suffix
  // (__upper -> tier 1, __bottom -> tier 0)
  void populate();

  // JYJ (2026-02-07) Load HB parasitic from tech definition
  // Reads hb_layer via R/C values from set_layer_rc -via hb_layer
  void loadHybridBondParasitics();

  // --- Instance tier queries ---

  // Get tier of a dbInst (returns cached value, -1 if unknown)
  int getInstTier(odb::dbInst* inst) const;

  // Get all instances on a given tier (0 or 1)
  const std::vector<odb::dbInst*>& getInstancesOnTier(int tier) const;

  // --- Buffer master mapping ---
  // JYJ (2026-02-06) Moved from HTreeBuilder::mapBufferMasterToTier()

  // Map a buffer master name to tier-specific variant.
  // First checks explicit bottom<->upper pair set via setTierBufferPair(),
  // then falls back to suffix substitution (_bottom <-> _upper).
  // e.g. with pair ("BUF_X4_bottom", "BUF_X4_upper"):
  //   "BUF_X4_bottom" + tier=1 -> "BUF_X4_upper"
  // e.g. suffix fallback (no pair set):
  //   "BUFx4_ASAP7_75t_R_bottom" + tier=1 -> "BUFx4_ASAP7_75t_R_upper"
  std::string getBufferForTier(const std::string& baseMaster,
                               int targetTier) const;

  // JYJ (2026-02-23) V32a: Register explicit bottom/upper buffer pair.
  // Bypasses suffix-only guessing when cell names differ across tech nodes.
  // Call after populate(), before tree building.
  // e.g. setTierBufferPair("BUF_X4_bottom", "BUF_X4_upper")
  void setTierBufferPair(const std::string& bottomBuf,
                         const std::string& upperBuf);

  // --- Dominant tier computation ---
  // JYJ (2026-02-06) Moved from HTreeBuilder::getDominantTierFrom*()

  // Majority voting over a set of ClockInst pointers
  int getDominantTier(const std::vector<ClockInst*>& insts) const;

  // Majority voting over locations using a location-to-sink map
  int getDominantTier(
      const std::vector<Point<double>>& locs,
      const std::map<Point<double>, ClockInst*>& locToSink) const;

  // Majority voting over all sinks in a Clock object
  int getDominantTierFromClock(const Clock& clock) const;

  // --- ClockInst tier management ---
  // JYJ (2026-02-07) Set tier on CTS-internal ClockInst through SSOT
  // (for buffers created during tree building, before writeDataToDb)
  void setClockInstTier(ClockInst& inst, int tier);

  // --- Cross-tier net detection ---

  // Returns true if a net connects instances on different tiers
  bool isCrossTierNet(odb::dbNet* net) const;

  // --- Per-tier wire parasitic ---
  // JYJ (2026-02-06) Extension point for 3D-aware parasitic estimation

  double getResPerDBU(int tier) const;
  double getCapPerDBU(int tier) const;
  void setWireRC(int tier, double resPerDBU, double capPerDBU);

  // --- Hybrid Bond (HBT) parasitic ---

  double getHbtResistance() const { return hbtRes_; }
  double getHbtCapacitance() const { return hbtCap_; }
  void setHbtParasitic(double resistance, double capacitance);

  // JYJ (2026-02-07) Convert HB R×C delay to equivalent wire distance
  // For use in clustering: penalize cross-tier pairs by this distance
  // Physics: HB lumped RC vs wire distributed RC (delay ∝ L²)
  // Formula: L = sqrt((R_HB × C_HB) / (r × c))
  // Returns: equivalent distance in μm (same units as wireResPerUnit)
  double getHbtEquivalentDistance(double wireResPerUnit,
                                  double wireCapPerUnit) const;

  // --- Pre-CTS skew targets (SG-CTS) ---
  // JYJ (2026-02-21) Load per-FF arrival time targets from LP solver
  void loadSkewTargets(const std::string& csv_path);
  double getSkewTarget(const std::string& ff_name) const;
  double getSkewTarget(odb::dbInst* inst) const;
  bool hasSkewTargets() const { return !skewTargetMap_.empty(); }
  int getSkewTargetCount() const
  {
    return static_cast<int>(skewTargetMap_.size());
  }

  // --- Statistics ---

  int getNumInstancesOnTier(int tier) const;
  void reportStats() const;

 private:
  // Determine tier from master name suffix
  static int tierFromMasterName(const std::string& name);

  // String suffix helper (moved from HTreeBuilder anonymous namespace)
  static bool hasSuffix(const std::string& name, const std::string& suffix);

  odb::dbDatabase* db_;
  utl::Logger* logger_;

  // Cached tier assignments (populated once from ODB)
  std::unordered_map<odb::dbInst*, int> instTierMap_;
  std::array<std::vector<odb::dbInst*>, 2> instsPerTier_;
  static const std::vector<odb::dbInst*> emptyInstVec_;

  // Per-tier wire RC (initially same for both tiers, can be split later)
  std::array<double, 2> resPerDBU_ = {0.0, 0.0};
  std::array<double, 2> capPerDBU_ = {0.0, 0.0};

  // JYJ (2026-02-23) V32a: Explicit bottom/upper buffer pair for cross-tech mapping.
  // Set via setTierBufferPair(); empty strings mean "use suffix fallback only".
  std::string bottomBufName_;
  std::string upperBufName_;

  // Hybrid Bond parasitic values (not TSV - this is face-to-face bonding)
  double hbtRes_ = 0.0;   // ohms
  double hbtCap_ = 0.0;   // farads

  // Pre-CTS skew targets: ff_instance_name -> target arrival offset (ns)
  std::unordered_map<std::string, double> skewTargetMap_;
};

}  // namespace cts
