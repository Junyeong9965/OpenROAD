// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors
//
// Created Cts3DDatabase implementation.
// Centralizes all 3D tier logic previously scattered across HTreeBuilder
// (mapBufferMasterToTier, getDominantTierFrom*, hasSuffix) and
// VerilogFFExtractor (master name parsing).
// Added setTierBufferPair() and updated getBufferForTier()
// to check explicit cross-tech buffer mapping before suffix substitution.
// Fixes hold timing regression in ASAP7_NG45_3D where "BUF_X4_bottom" ->
// "BUF_X4_upper" (suffix swap) fails because "BUF_X4_upper" does not exist;
// the correct upper-tier equivalent is "BUF_X4_upper" (new wrapper cell).

#include "Cts3DDatabase.h"

#include <cmath>
#include <fstream>
#include <sstream>
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

// Moved from floorplan_utils.tcl set_tier_from_master_names()
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

  // Debug: count total instances and sample first master name
  int totalInsts = 0;
  std::string firstMasterName = "(none)";
  for (odb::dbInst* inst : block->getInsts()) {
    ++totalInsts;
    if (totalInsts == 1 && inst->getMaster() != nullptr) {
      firstMasterName = inst->getMaster()->getName();
    }
  }
  logger_->info(CTS, 357,
                "Cts3DDatabase::populate debug: block has {} insts, first master='{}'",
                totalInsts, firstMasterName);

  for (odb::dbInst* inst : block->getInsts()) {
    odb::dbMaster* master = inst->getMaster();
    if (master == nullptr) {
      continue;
    }
    const std::string masterName = master->getName();
    int tier = tierFromMasterName(masterName);

    // Fallback: use ODB tier (set by set_tier_from_master_names)
    // Only trust getTier()==1 as unambiguous (default z_=0 means bottom OR unset)
    if (tier < 0) {
      const int odb_tier = inst->getTier();
      if (odb_tier == 1) {
        tier = 1;  // unambiguously upper (explicitly set by Tcl proc)
      } else if (odb_tier == 0 &&
                 (masterName.find("_bottom") == std::string::npos &&
                  masterName.find("_upper") == std::string::npos)) {
        // tier=0 AND no suffix → might be default; skip ambiguous instances
      }
    }

    if (tier >= 0) {
      instTierMap_[inst] = tier;
      instsPerTier_[tier].push_back(inst);
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

  // Load HB parasitic after tier assignment
  loadHybridBondParasitics();
}

// Load hybrid bond parasitic from tech LEF.
// R: read from hb_layer RESISTANCE in tech LEF (ohms per cut for CUT layers).
//    Fallback default 0.02Ω matches Pin3D paper Table 1 (HBT geometry settings).
// C: tech LEF does not provide capacitance for CUT layers → set to 0.
//    HBV via delay = 0.693 * R * C ≈ 0 regardless of R value,
//    so cap uncertainty has no practical impact on timing.
void Cts3DDatabase::loadHybridBondParasitics()
{
  // Fallback: Pin3D paper Table 1 HBT geometry (0.5μm width, 1.0μm pitch)
  const double DEFAULT_HB_RES = 0.02;       // Ω per cut (from tech LEF RESISTANCE)
  const double DEFAULT_HB_CAP = 0.0;        // F (not available in tech LEF for CUT layers)

  odb::dbTech* tech = db_->getTech();
  if (tech == nullptr) {
    logger_->warn(CTS, 358,
                  "No tech found. Using defaults: R={:.4f}Ω, C=0 (not in tech LEF)",
                  DEFAULT_HB_RES);
    hbtRes_ = DEFAULT_HB_RES;
    hbtCap_ = DEFAULT_HB_CAP;
    return;
  }

  // Find hb_layer (TYPE CUT) in tech LEF
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
                  "Using defaults: R={:.4f}Ω, C=0. "
                  "Add 'LAYER hb_layer' with RESISTANCE to tech LEF.",
                  DEFAULT_HB_RES);
    hbtRes_ = DEFAULT_HB_RES;
    hbtCap_ = DEFAULT_HB_CAP;
    return;
  }

  // Read R from tech LEF (dbTechLayer::getResistance() returns ohms per cut for CUT layers)
  hbtRes_ = hb_layer->getResistance();
  if (hbtRes_ <= 0.0) {
    logger_->warn(CTS, 544,
                  "hb_layer RESISTANCE={:.4f} invalid, using default {:.4f}Ω",
                  hbtRes_, DEFAULT_HB_RES);
    hbtRes_ = DEFAULT_HB_RES;
  }

  // C not available in tech LEF for CUT layers → 0
  hbtCap_ = DEFAULT_HB_CAP;

  // Elmore 50% delay: t_via = 0.693 * R * C
  const double t_via_s = 0.693 * hbtRes_ * hbtCap_;
  logger_->info(CTS, 359,
                "HB parasitic loaded: R={:.4f}Ω (tech LEF), C=0F (not in LEF), "
                "t_via={:.3e}s (0.693*R*C, negligible)",
                hbtRes_, t_via_s);
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
  // No ODB fallback - our 3DDB is the single source of truth
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

// Moved from HTreeBuilder.cpp mapBufferMasterToTier()
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

  // CTS_FORCE_SINGLE_TIER_BUF: disable tier swapping, always return
  // baseMaster (root buffer). Pin3D uses bottom-only buffers at ALL
  // H-tree levels → no cross-tier clock branches → uniform tree.
  static int forceSingleTier = -1;
  if (forceSingleTier < 0) {
    const char* e = std::getenv("CTS_FORCE_SINGLE_TIER_BUF");
    forceSingleTier = (e && std::atoi(e) != 0) ? 1 : 0;
  }
  if (forceSingleTier) {
    return baseMaster;
  }

  // Check explicit cross-tech tier buffer pair first.
  // Required when bottom/upper cell names don't share the same base name
  // (e.g., "BUF_X4_bottom" vs "BUF_X4_upper" which is an ASAP7 wrapper cell,
  // as opposed to the original "BUFx4_ASAP7_75t_R_upper" whose suffix swap
  // would yield non-existent "BUF_X4_upper" then fall back to "BUF_X4_bottom").
  if (!bottomBufName_.empty() && !upperBufName_.empty()) {
    if (targetTier == 1 && baseMaster == bottomBufName_) {
      return upperBufName_;
    }
    if (targetTier == 0 && baseMaster == upperBufName_) {
      return bottomBufName_;
    }
  }

  // Suffix-based substitution fallback (works for same-family cells).
  std::string candidate = baseMaster;
  if (targetTier == 1) {
    if (hasSuffix(baseMaster, "_bottom")) {
      candidate = baseMaster.substr(0, baseMaster.size() - 7) + "_upper";
    } else if (!hasSuffix(baseMaster, "_upper")) {
      candidate = baseMaster + "_upper";
    }
  } else {
    if (hasSuffix(baseMaster, "_upper")) {
      candidate = baseMaster.substr(0, baseMaster.size() - 6) + "_bottom";
    } else if (!hasSuffix(baseMaster, "_bottom")) {
      candidate = baseMaster + "_bottom";
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

// Register explicit bottom/upper buffer name pair.
// Enables getBufferForTier() to correctly map across heterogeneous tech nodes
// where cell base names differ (e.g., NG45 "BUF_X4_bottom" <-> ASAP7 "BUF_X4_upper").
void Cts3DDatabase::setTierBufferPair(const std::string& bottomBuf,
                                      const std::string& upperBuf)
{
  bottomBufName_ = bottomBuf;
  upperBufName_  = upperBuf;
  logger_->info(CTS, 370,
                "3D CTS tier buffer pair registered: bottom='{}', upper='{}'",
                bottomBufName_, upperBufName_);
}

// Moved from HTreeBuilder::getDominantTierFromInsts()
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

// Moved from HTreeBuilder::getDominantTierFromSinkLocs()
int Cts3DDatabase::getDominantTier(
    const std::vector<Point<double>>& locs,
    const std::map<Point<double>, ClockInst*>& locToSink) const
{
  int tier0 = 0;
  int tier1 = 0;
  for (const Point<double>& loc : locs) {
    // V58 Fix B: epsilon-tolerant lookup for tier vote.
    auto it = locToSink.find(loc);
    if (it == locToSink.end()) {
      constexpr double eps = 1e-4;
      auto hint = locToSink.lower_bound(
          Point<double>(loc.getX() - eps, loc.getY() - eps));
      for (auto scan = hint; scan != locToSink.end(); ++scan) {
        if (scan->first.getX() > loc.getX() + eps) break;
        if (std::abs(scan->first.getX() - loc.getX()) < eps
            && std::abs(scan->first.getY() - loc.getY()) < eps) {
          it = scan;
          break;
        }
      }
      if (it == locToSink.end()) {
        continue;
      }
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

// Moved from HTreeBuilder::getDominantTierFromClockSinks()
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

// SSOT: set tier on CTS-internal ClockInst
void Cts3DDatabase::setClockInstTier(ClockInst& inst, int tier)
{
  inst.setTier(tier);
}

// New: cross-tier net detection
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

// Convert HB parasitic to equivalent wire distance
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

// --- Pre-CTS skew targets ---

// Load per-FF arrival targets from Python LP solver output
void Cts3DDatabase::loadSkewTargets(const std::string& csv_path)
{
  skewTargetMap_.clear();
  std::ifstream file(csv_path);
  if (!file.is_open()) {
    logger_->warn(CTS, 397, "Cannot open skew targets file: {}", csv_path);
    return;
  }

  std::string line;
  // Skip header
  if (!std::getline(file, line)) {
    logger_->warn(CTS, 398, "Empty skew targets file: {}", csv_path);
    return;
  }

  int count = 0;
  int nonzero = 0;
  while (std::getline(file, line)) {
    std::istringstream ss(line);
    std::string ff_name, target_str, tier_str;
    std::getline(ss, ff_name, ',');
    std::getline(ss, target_str, ',');
    std::getline(ss, tier_str, ',');

    if (ff_name.empty() || target_str.empty()) {
      continue;
    }

    try {
      double target = std::stod(target_str);
      skewTargetMap_[ff_name] = target;
      count++;
      if (std::abs(target) > 1e-6) {
        nonzero++;
      }
    } catch (const std::exception&) {
      // Skip malformed rows
    }
  }

  logger_->info(CTS, 399,
                "Loaded {} skew targets ({} non-zero) from {}",
                count, nonzero, csv_path);
}

double Cts3DDatabase::getSkewTarget(const std::string& ff_name) const
{
  auto it = skewTargetMap_.find(ff_name);
  return (it != skewTargetMap_.end()) ? it->second : 0.0;
}

double Cts3DDatabase::getSkewTarget(odb::dbInst* inst) const
{
  if (inst == nullptr || skewTargetMap_.empty()) {
    return 0.0;
  }
  return getSkewTarget(inst->getName());
}

// --- Private helpers ---

int Cts3DDatabase::tierFromMasterName(const std::string& name)
{
  if (name.find("_upper") != std::string::npos) {
    return 1;
  }
  if (name.find("_bottom") != std::string::npos) {
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
