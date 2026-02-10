// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors
//
// JYJ (2026-02-06) Created Cts3DDatabase: centralized 3D tier management
// for True 3D CTS. Replaces scattered tier logic in HTreeBuilder and
// VerilogFFExtractor with a single source of truth.

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

  // Map a buffer master name to tier-specific variant
  // e.g. "BUFx4_ASAP7_75t_R__bottom" + tier=1 -> "BUFx4_ASAP7_75t_R__upper"
  std::string getBufferForTier(const std::string& baseMaster,
                               int targetTier) const;

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

  // Hybrid Bond parasitic values (not TSV - this is face-to-face bonding)
  double hbtRes_ = 0.0;   // ohms
  double hbtCap_ = 0.0;   // farads
};

}  // namespace cts
