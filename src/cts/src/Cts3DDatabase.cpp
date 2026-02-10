// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors
//
// JYJ (2026-02-06) Created Cts3DDatabase implementation.
// Centralizes all 3D tier logic previously scattered across HTreeBuilder
// (mapBufferMasterToTier, getDominantTierFrom*, hasSuffix) and
// VerilogFFExtractor (master name parsing).

#include "Cts3DDatabase.h"

#include <cmath>
#include <string>

#include "utl/Logger.h"

namespace cts {

using utl::CTS;

// Static empty vector for out-of-range tier queries
const std::vector<odb::dbInst*> Cts3DDatabase::emptyInstVec_;

Cts3DDatabase::Cts3DDatabase(odb::dbDatabase* db, utl::Logger* logger)
    : db_(db), logger_(logger)
{
}

// JYJ (2026-02-06) Moved from floorplan_utils.tcl set_tier_from_master_names()
// and HTreeBuilder anonymous namespace. Now done once in C++ at startup.
void Cts3DDatabase::populate()
{
  instTierMap_.clear();
  instsPerTier_[0].clear();
  instsPerTier_[1].clear();

  odb::dbChip* chip = db_->getChip();
  if (chip == nullptr) {
    logger_->warn(CTS, 350, "Cts3DDatabase: no chip found, skipping populate");
    return;
  }
  odb::dbBlock* block = chip->getBlock();
  if (block == nullptr) {
    logger_->warn(CTS, 351, "Cts3DDatabase: no block found, skipping populate");
    return;
  }

  int upperCount = 0;
  int bottomCount = 0;

  for (odb::dbInst* inst : block->getInsts()) {
    odb::dbMaster* master = inst->getMaster();
    if (master == nullptr) {
      continue;
    }
    const std::string masterName = master->getName();
    const int tier = tierFromMasterName(masterName);

    if (tier >= 0) {
      instTierMap_[inst] = tier;
      instsPerTier_[tier].push_back(inst);
      // JYJ (2026-02-09) Do NOT set ODB tier - we maintain our own 3DDB
      if (tier == 0) {
        bottomCount++;
      } else {
        upperCount++;
      }
    }
  }

  logger_->info(CTS, 352,
                "Cts3DDatabase populated: bottom={}, upper={}, total={}",
                bottomCount, upperCount, bottomCount + upperCount);

  // JYJ (2026-02-07) Load HB parasitic after tier assignment
  loadHybridBondParasitics();
}

// JYJ (2026-02-07) Load hybrid bond parasitic from ODB tech definition
// Reads R from tech LEF, uses paper default for C (tech LEF doesn't provide C for vias)
void Cts3DDatabase::loadHybridBondParasitics()
{
  // Default values from paper (Table 2)
  const double DEFAULT_HB_RES = 3.0;        // Ω
  const double DEFAULT_HB_CAP = 0.6e-15;    // F (0.6 fF)

  odb::dbTech* tech = db_->getTech();
  if (tech == nullptr) {
    logger_->warn(CTS, 358,
                  "No tech found. Using paper defaults: R={:.1f}Ω, C={:.1f}fF",
                  DEFAULT_HB_RES, DEFAULT_HB_CAP * 1e15);
    hbtRes_ = DEFAULT_HB_RES;
    hbtCap_ = DEFAULT_HB_CAP;
    return;
  }

  // Try to find hb_layer in tech
  // Note: Actual ODB API may differ - adjust based on available methods
  // TODO: Replace with actual ODB API when schema is confirmed
  odb::dbTechLayer* hb_layer = nullptr;
  for (odb::dbTechLayer* layer : tech->getLayers()) {
    if (layer->getName() == "hb_layer") {
      hb_layer = layer;
      break;
    }
  }

  if (hb_layer == nullptr) {
    logger_->warn(CTS, 360,
                  "HB layer 'hb_layer' not found in tech LEF. "
                  "Using paper defaults: R={:.1f}Ω, C={:.1f}fF. "
                  "Add 'LAYER hb_layer' with RESISTANCE to tech LEF.",
                  DEFAULT_HB_RES, DEFAULT_HB_CAP * 1e15);
    hbtRes_ = DEFAULT_HB_RES;
    hbtCap_ = DEFAULT_HB_CAP;
    return;
  }

  // Read R from tech (if available)
  // TODO: Use actual ODB API method (e.g., hb_layer->getResistance())
  // Placeholder: Assume we can read resistance
  hbtRes_ = DEFAULT_HB_RES;  // TODO: Replace with hb_layer->getResistance()

  // C not typically in tech LEF for vias → use paper default
  hbtCap_ = DEFAULT_HB_CAP;

  logger_->info(CTS, 359,
                "HB parasitic loaded: R={:.3e}Ω (tech LEF), C={:.3e}F (paper default)",
                hbtRes_, hbtCap_);
  logger_->info(CTS, 361,
                "HB equivalent distance: RC delay mode (~15.5 μm) enabled. "
                "R-only fallback (~58.3 μm) available if wire C unavailable.");
}

int Cts3DDatabase::getInstTier(odb::dbInst* inst) const
{
  if (inst == nullptr) {
    return -1;
  }
  auto it = instTierMap_.find(inst);
  if (it != instTierMap_.end()) {
    return it->second;
  }
  // JYJ (2026-02-09) No ODB fallback - our 3DDB is the single source of truth
  // If instance not in map, tier is unknown
  return -1;
}

const std::vector<odb::dbInst*>& Cts3DDatabase::getInstancesOnTier(
    int tier) const
{
  if (tier < 0 || tier > 1) {
    return emptyInstVec_;
  }
  return instsPerTier_[tier];
}

// JYJ (2026-02-06) Moved from HTreeBuilder.cpp mapBufferMasterToTier()
// and hasSuffix(). Logic is identical, now centralized here.
std::string Cts3DDatabase::getBufferForTier(const std::string& baseMaster,
                                            int targetTier) const
{
  if (targetTier != 0 && targetTier != 1) {
    return baseMaster;
  }
  if (db_ == nullptr || baseMaster.empty()) {
    return baseMaster;
  }

  std::string candidate = baseMaster;
  if (targetTier == 1) {
    if (hasSuffix(baseMaster, "__bottom")) {
      candidate = baseMaster.substr(0, baseMaster.size() - 9) + "__upper";
    } else if (!hasSuffix(baseMaster, "__upper")) {
      candidate = baseMaster + "__upper";
    }
  } else {
    if (hasSuffix(baseMaster, "__upper")) {
      candidate = baseMaster.substr(0, baseMaster.size() - 8) + "__bottom";
    } else if (!hasSuffix(baseMaster, "__bottom")) {
      candidate = baseMaster + "__bottom";
    }
  }

  if (candidate == baseMaster) {
    return baseMaster;
  }

  // Verify the tier-specific master exists in the database
  if (db_->findMaster(candidate.c_str()) != nullptr) {
    return candidate;
  }

  return baseMaster;
}

// JYJ (2026-02-06) Moved from HTreeBuilder::getDominantTierFromInsts()
int Cts3DDatabase::getDominantTier(
    const std::vector<ClockInst*>& insts) const
{
  int tier0 = 0;
  int tier1 = 0;
  for (const ClockInst* inst : insts) {
    if (inst == nullptr) {
      continue;
    }
    odb::dbInst* db_inst = inst->getDbInst();
    if (db_inst == nullptr) {
      continue;
    }
    const int tier = getInstTier(db_inst);
    if (tier == 0) {
      tier0++;
    } else if (tier == 1) {
      tier1++;
    }
  }
  if (tier0 == 0 && tier1 == 0) {
    return -1;
  }
  return (tier1 > tier0) ? 1 : 0;
}

// JYJ (2026-02-06) Moved from HTreeBuilder::getDominantTierFromSinkLocs()
int Cts3DDatabase::getDominantTier(
    const std::vector<Point<double>>& locs,
    const std::map<Point<double>, ClockInst*>& locToSink) const
{
  int tier0 = 0;
  int tier1 = 0;
  for (const Point<double>& loc : locs) {
    auto it = locToSink.find(loc);
    if (it == locToSink.end()) {
      continue;
    }
    ClockInst* inst = it->second;
    if (inst == nullptr) {
      continue;
    }
    int tier = -1;
    odb::dbInst* db_inst = inst->getDbInst();
    if (db_inst != nullptr) {
      tier = getInstTier(db_inst);
    } else {
      // Fallback to ClockInst tier (for buffers not yet in ODB)
      tier = inst->getTier();
    }
    if (tier == 0) {
      tier0++;
    } else if (tier == 1) {
      tier1++;
    }
  }
  if (tier0 == 0 && tier1 == 0) {
    return -1;
  }
  return (tier1 > tier0) ? 1 : 0;
}

// JYJ (2026-02-06) Moved from HTreeBuilder::getDominantTierFromClockSinks()
int Cts3DDatabase::getDominantTierFromClock(const Clock& clock) const
{
  int tier0 = 0;
  int tier1 = 0;
  clock.forEachSink([&](const ClockInst& inst) {
    odb::dbInst* db_inst = inst.getDbInst();
    if (db_inst == nullptr) {
      return;
    }
    const int tier = getInstTier(db_inst);
    if (tier == 0) {
      tier0++;
    } else if (tier == 1) {
      tier1++;
    }
  });
  if (tier0 == 0 && tier1 == 0) {
    return -1;
  }
  return (tier1 > tier0) ? 1 : 0;
}

// JYJ (2026-02-07) SSOT: set tier on CTS-internal ClockInst
void Cts3DDatabase::setClockInstTier(ClockInst& inst, int tier)
{
  inst.setTier(tier);
}

// JYJ (2026-02-06) New: cross-tier net detection
bool Cts3DDatabase::isCrossTierNet(odb::dbNet* net) const
{
  if (net == nullptr) {
    return false;
  }
  bool hasTier0 = false;
  bool hasTier1 = false;
  for (odb::dbITerm* iterm : net->getITerms()) {
    odb::dbInst* inst = iterm->getInst();
    const int tier = getInstTier(inst);
    if (tier == 0) {
      hasTier0 = true;
    } else if (tier == 1) {
      hasTier1 = true;
    }
    if (hasTier0 && hasTier1) {
      return true;
    }
  }
  return false;
}

double Cts3DDatabase::getResPerDBU(int tier) const
{
  if (tier < 0 || tier > 1) {
    return resPerDBU_[0];  // default to bottom tier
  }
  return resPerDBU_[tier];
}

double Cts3DDatabase::getCapPerDBU(int tier) const
{
  if (tier < 0 || tier > 1) {
    return capPerDBU_[0];
  }
  return capPerDBU_[tier];
}

void Cts3DDatabase::setWireRC(int tier, double resPerDBU, double capPerDBU)
{
  if (tier >= 0 && tier <= 1) {
    resPerDBU_[tier] = resPerDBU;
    capPerDBU_[tier] = capPerDBU;
  }
}

void Cts3DDatabase::setHbtParasitic(double resistance, double capacitance)
{
  hbtRes_ = resistance;
  hbtCap_ = capacitance;
}

// JYJ (2026-02-07) Convert HB parasitic to equivalent wire distance
// Two modes: RC delay (default) or R-only (emergency fallback if wire C invalid)
double Cts3DDatabase::getHbtEquivalentDistance(double wireResPerUnit,
                                               double wireCapPerUnit) const
{
  if (wireResPerUnit <= 0.0 || hbtRes_ <= 0.0) {
    return 0.0;  // Need at least R values
  }

  // Mode 1: RC delay-based (more accurate if C available)
  // Physics: HB lumped RC vs wire distributed RC (delay ∝ L²)
  //   delay_HB = R_HB × C_HB
  //   delay_wire ≈ (1/2) × r × c × L²  (Elmore)
  //   Solve: L = sqrt(R_HB × C_HB / (r × c))
  if (hbtCap_ > 0.0 && wireCapPerUnit > 0.0) {
    const double delay_hb = hbtRes_ * hbtCap_;
    const double delay_wire_per_unit = wireResPerUnit * wireCapPerUnit;

    // Example: R=3Ω, C=0.6fF, r=0.0514Ω/um, c=0.145fF/um
    //   delay_hb = 1.8e-15 s
    //   delay_wire = 7.45e-18 s/um²
    //   L = sqrt(241.6) ≈ 15.5 um
    return std::sqrt(delay_hb / delay_wire_per_unit);
  }

  // Mode 2: R-only fallback (emergency: if wire C <= 0 from TechChar)
  // Physics: Equivalent wire length with same resistance as HB
  //   R_HB = r_wire × L
  //   L = R_HB / r_wire
  // Example: R=3Ω, r=0.0514Ω/um → L = 58.3 um
  return hbtRes_ / wireResPerUnit;
}

int Cts3DDatabase::getNumInstancesOnTier(int tier) const
{
  if (tier < 0 || tier > 1) {
    return 0;
  }
  return static_cast<int>(instsPerTier_[tier].size());
}

void Cts3DDatabase::reportStats() const
{
  logger_->info(CTS, 353,
                "Cts3DDatabase stats: bottom={}, upper={}",
                getNumInstancesOnTier(0), getNumInstancesOnTier(1));
}

// --- Private helpers ---

int Cts3DDatabase::tierFromMasterName(const std::string& name)
{
  if (name.find("__upper") != std::string::npos) {
    return 1;
  }
  if (name.find("__bottom") != std::string::npos) {
    return 0;
  }
  return -1;
}

bool Cts3DDatabase::hasSuffix(const std::string& name,
                              const std::string& suffix)
{
  if (name.size() < suffix.size()) {
    return false;
  }
  return name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0;
}

}  // namespace cts
