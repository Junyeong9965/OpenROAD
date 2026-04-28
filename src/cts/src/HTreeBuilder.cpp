// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "HTreeBuilder.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "Clustering.h"
#include "Cts3DDatabase.h"  // Added for 3D tier management
#include "SinkClustering.h"
#include "TarpClustering.h"  // TARP clustering
#include "TechChar.h"
#include "TreeBuilder.h"
#include "Util.h"
#include "odb/db.h"
#include "odb/isotropy.h"
#include "utl/Logger.h"

namespace cts {

using utl::CTS;

namespace {

bool isNamedHardMacroMaster(const odb::dbMaster* master)
{
  if (master == nullptr) {
    return false;
  }
  std::string name = master->getName();
  std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) {
    return std::tolower(c);
  });
  return name.rfind("fakeram", 0) == 0 || name.rfind("fakeregfile", 0) == 0
         || name.rfind("sram_", 0) == 0;
}

bool isMacroBlockInst(const odb::dbInst* inst)
{
  if (inst == nullptr) {
    return false;
  }
  if (inst->isBlock()) {
    return true;
  }
  const odb::dbMaster* master = inst->getMaster();
  return (master != nullptr && master->isBlock())
         || isNamedHardMacroMaster(master);
}

}  // namespace

// Removed anonymous namespace containing hasSuffix() and
// mapBufferMasterToTier() — moved to Cts3DDatabase::hasSuffix() and
// Cts3DDatabase::getBufferForTier()

Point<double> HTreeBuilder::legalizeOneBuffer(Point<double> bufferLoc,
                                              const std::string& bufferName)
{
  Point<double> legalLoc
      = TreeBuilder::legalizeOneBuffer(bufferLoc, bufferName);
  return resolveLocationCollision(legalLoc);
}

// Step 3: Skew-aware clustering via skew target penalty
// HB penalty removed — vias (including HB) don't get distance penalty in 2D CTS
// either, and STA already accounts for HB RC in timing. Clustering should remain
// purely spatial + skew-target-driven.
// Removed skewTargetBeta_ dead code block.
// computeDist() is called with candidate branching points (not in
// mapLocationToSink_) → find() always failed → beta penalty never fired.
// Target-aware clustering works via CTS_CKMEANS_TARGET_BETA in
// refineBranchingPointsWithClustering() → Clustering::calcDist().
double HTreeBuilder::computeDist(const Point<double>& x, const Point<double>& y)
{
  return TreeBuilder::computeDist(x, y);
}

// V58 Fix B: Epsilon-tolerant lookup for mapLocationToSink_.
// Falls back to linear scan within eps when exact find() fails.
// Matches the epsilon pattern already used in leaf-sink connection code.
ClockInst* HTreeBuilder::findSinkEps(const Point<double>& loc,
                                      double eps) const
{
  auto it = mapLocationToSink_.find(loc);
  if (it != mapLocationToSink_.end()) {
    return it->second;
  }
  auto hint = mapLocationToSink_.lower_bound(
      Point<double>(loc.getX() - eps, loc.getY() - eps));
  for (auto scan = hint; scan != mapLocationToSink_.end(); ++scan) {
    if (scan->first.getX() > loc.getX() + eps) break;
    if (std::abs(scan->first.getX() - loc.getX()) < eps
        && std::abs(scan->first.getY() - loc.getY()) < eps) {
      return scan->second;
    }
  }
  return nullptr;
}

void HTreeBuilder::preSinkClustering(
    const std::vector<std::pair<float, float>>& sinks,
    const std::vector<const ClockInst*>& sinkInsts,
    const float maxDiameter,
    const unsigned clusterSize,
    const bool secondLevel)
{
  bool maxDiameterSet = (type_ == TreeType::MacroTree)
                            ? options_->isMacroMaxDiameterSet()
                            : options_->isMaxDiameterSet();
  bool clusterSizeSet = (type_ == TreeType::MacroTree)
                            ? options_->isMacroSinkClusteringSizeSet()
                            : options_->isSinkClusteringSizeSet();

  unsigned min_clustering_sinks = (type_ == TreeType::MacroTree)
                                      ? min_clustering_macro_sinks_
                                      : min_clustering_sinks_;

  const std::vector<std::pair<float, float>>& points = sinks;
  if (!secondLevel) {
    clock_.forEachSink([&](ClockInst& inst) {
      // FloatFix: (float) → (double) to prevent map key mismatch → CTS-0080
      const Point<double> normLocation((double) inst.getX() / wireSegmentUnit_,
                                       (double) inst.getY() / wireSegmentUnit_);
      mapLocationToSink_[normLocation] = &inst;
      if (!fuzzyEqual(inst.getInsertionDelay(), 0.0, 1e-6)) {
        setSinkInsertionDelay(normLocation,
                              inst.getInsertionDelay() / wireSegmentUnit_);
        // clang-format off
	debugPrint(logger_, CTS, "clustering", 1, "sink {} has insDelay {} at {}",
		   inst.getName(), getSinkInsertionDelay(normLocation),
		   normLocation);
        // clang-format on
      }
    });
  }

  if (sinks.size() <= min_clustering_sinks
      || !(options_->getSinkClustering())) {
    topLevelSinksClustered_ = sinks;
    return;
  }

  // Hilbert-calibrated TARP.
  // Run Hilbert grid search as lightweight pre-pass to find optimal
  // cluster granularity (bestGroupSize + bestDiameter), then pass
  // to TARP FM. This ensures TARP operates at the right spatial
  // granularity for each design without per-design manual tuning.
  if (enableTarp_ && !secondLevel
      && cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()
      && sinks.size() >= 50u) {

    // Step 1: Lightweight Hilbert grid search (no leaf buffer creation)
    SinkClustering calibration(options_, techChar_, this);
    const unsigned numPts = points.size();
    for (unsigned i = 0; i < numPts; ++i) {
      calibration.addPoint(points[i].first, points[i].second);
      if (sinkInsts[i]->getInputCap() == 0) {
        calibration.addCap(options_->getSinkBufferInputCap());
      } else {
        calibration.addCap(sinkInsts[i]->getInputCap());
      }
    }

    unsigned bestGS = 0;
    float bestDiam = 0.0f;
    // Grid search: groupSize ∈ {10,20,30} × diameter ∈ {50,100,200}μm
    for (unsigned diam : clusterDiameters()) {
      for (unsigned gs : options_->getSinkClusteringSizes()) {
        float maxDiam = diam * static_cast<float>(options_->getDbUnits())
                       / wireSegmentUnit_;
        calibration.run(gs, maxDiam, wireSegmentUnit_, bestGS, bestDiam);
      }
    }

    logger_->info(CTS, 845,
        "TARP Hilbert calibration: bestGroupSize={}, bestDiameter={:.0f}um",
        bestGS, bestDiam / options_->getDbUnits() * wireSegmentUnit_);

    // Step 2: Pass calibrated params to TARP (env var overrides take priority)
    runTarpClustering(points, sinkInsts, bestGS, bestDiam);
    return;
  }

  SinkClustering matching(options_, techChar_, this);
  const unsigned numPoints = points.size();

  for (int pointIdx = 0; pointIdx < numPoints; ++pointIdx) {
    const std::pair<float, float>& point = points[pointIdx];
    matching.addPoint(point.first, point.second);
    if (sinkInsts[pointIdx]->getInputCap() == 0) {
      // Comes here in second level since first level buf cap is not set
      matching.addCap(options_->getSinkBufferInputCap());
    } else {
      matching.addCap(sinkInsts[pointIdx]->getInputCap());
    }
  }

  unsigned bestClusterSize = 0;
  float bestDiameter = 0.0;
  if (clusterSizeSet && maxDiameterSet) {
    // clang-format off
      debugPrint(logger_, CTS, "clustering", 1, "**** match.run({}, {}, {}) ****",
                 clusterSize, maxDiameter, wireSegmentUnit_);
    // clang-format on
    matching.run(clusterSize,
                 maxDiameter,
                 wireSegmentUnit_,
                 bestClusterSize,
                 bestDiameter);
  } else if (!clusterSizeSet && maxDiameterSet) {
    // only diameter is set, try clustering sizes of 10, 20 and 30
    for (unsigned clusterSize2 : options_->getSinkClusteringSizes()) {
      // clang-format off
      debugPrint(logger_, CTS, "clustering", 1, "**** match.run({}, {}, {}) ****",
                 clusterSize2, maxDiameter, wireSegmentUnit_);
      // clang-format on
      matching.run(clusterSize2,
                   maxDiameter,
                   wireSegmentUnit_,
                   bestClusterSize,
                   bestDiameter);
    }
  } else if (clusterSizeSet && !maxDiameterSet) {
    // only clustering size is set, try diameters of 50, 100 and 200 um
    for (unsigned clusterDiameter2 : options_->getSinkClusteringDiameters()) {
      // clang-format off
      debugPrint(logger_, CTS, "clustering", 1, "**** match.run({}, {}, {}) ****",
                 clusterSize, clusterDiameter2, wireSegmentUnit_);
      // clang-format on
      float maxDiameter2 = clusterDiameter2 * (float) options_->getDbUnits()
                           / wireSegmentUnit_;
      matching.run(clusterSize,
                   maxDiameter2,
                   wireSegmentUnit_,
                   bestClusterSize,
                   bestDiameter);
    }
  } else {  // neighther clustering size nor diameter is set
    // try diameters of 50, 100 and 200 um
    for (unsigned clusterDiameter2 : clusterDiameters()) {
      // try clustering sizes of 10, 20 and 30
      for (unsigned clusterSize2 : options_->getSinkClusteringSizes()) {
        // clang-format off
        debugPrint(logger_, CTS, "clustering", 1, "**** match.run({}, {}, {}) ****",
                   clusterSize2, clusterDiameter2, wireSegmentUnit_);
        // clang-format on
        float maxDiameter2 = clusterDiameter2 * (float) options_->getDbUnits()
                             / wireSegmentUnit_;
        matching.run(clusterSize2,
                     maxDiameter2,
                     wireSegmentUnit_,
                     bestClusterSize,
                     bestDiameter);
      }
    }
  }

  if (clusterSizeSet || maxDiameterSet) {
    logger_->info(
        CTS,
        204,
        "A clustering solution was found from clustering size of {} "
        "and clustering diameter of {:0.0f}.",
        bestClusterSize,
        std::round(bestDiameter / options_->getDbUnits() * wireSegmentUnit_));
    logger_->info(CTS,
                  205,
                  "Better solution may be possible if either "
                  "-sink_clustering_size, -sink_clustering_max_diameter, or "
                  "both options are omitted to enable automatic clustering.");
  } else {
    logger_->info(
        CTS,
        206,
        "Best clustering solution was found from clustering size of {} "
        "and clustering diameter of {:0.0f}.",
        bestClusterSize,
        std::round(bestDiameter / options_->getDbUnits() * wireSegmentUnit_));
  }

  unsigned clusterCount = 0;

  // Tier-aware clustering: detect cross-tier clusters
  // Calculate HB equivalent distance for cross-tier penalty awareness
  const double wire_res_per_unit = cts3dDb_->getResPerDBU(0) * wireSegmentUnit_;
  const double wire_cap_per_unit = cts3dDb_->getCapPerDBU(0) * wireSegmentUnit_;
  const double hb_equivalent_dist = cts3dDb_->getHbtEquivalentDistance(
      wire_res_per_unit, wire_cap_per_unit);

  // Note: Current k-means clustering uses geometric distance only.
  // Possible extension: inject HB equivalent distance into SinkClustering/CKMeans
  // as penalty when computing distance between sinks on different tiers.
  // This would actively discourage cross-tier clusters during k-means.

  std::vector<std::pair<float, float>> newSinkLocations;

  // Helper lambda to create one leaf buffer for a group of FFs.
  // This is factored out so that sub-clustering can call it multiple times per cluster.
  int splitCount_2way = 0, splitCount_3way = 0;  // sub-cluster stats
  auto createLeafBuffer = [&](const std::vector<ClockInst*>& subInsts,
                               const std::vector<unsigned>& subPointIdxs,
                               unsigned& leafIdx) {
    if (subInsts.empty()) return;
    float xSum = 0, ySum = 0;
    for (unsigned pi : subPointIdxs) {
      xSum += points[pi].first;
      ySum += points[pi].second;
    }
    const float normCenterX = xSum / subPointIdxs.size();
    const float normCenterY = ySum / subPointIdxs.size();
    Point<double> center((double) normCenterX, (double) normCenterY);

    const int target_tier = cts3dDb_->getDominantTier(subInsts);
    const std::string sink_buffer
        = cts3dDb_->getBufferForTier(options_->getSinkBuffer(), target_tier);

    const std::string baseName = tierPrefix_ + (secondLevel ? "clkbuf_leaf2_" : "clkbuf_leaf_");
    Point<double> rootBufLoc = legalizeOneBuffer(center, sink_buffer);
    commitMoveLoc(center, rootBufLoc);

    ClockInst& rootBuffer = clock_.addClockBuffer(
        baseName + std::to_string(leafIdx),
        sink_buffer,
        rootBufLoc.getX() * wireSegmentUnit_,
        rootBufLoc.getY() * wireSegmentUnit_);
    cts3dDb_->setClockInstTier(rootBuffer, target_tier);

    if (!secondLevel) {
      addFirstLevelSinkDriver(&rootBuffer);
    } else {
      addSecondLevelSinkDriver(&rootBuffer);
    }

    const std::string netBaseName = tierPrefix_ + (secondLevel ? "clknet_leaf2_" : "clknet_leaf_");
    ClockSubNet& clockSubNet
        = clock_.addSubNet(netBaseName + std::to_string(leafIdx));
    clockSubNet.addInst(rootBuffer);
    for (ClockInst* ci : subInsts) {
      clockSubNet.addInst(*ci);
    }
    if (!secondLevel) {
      clockSubNet.setLeafLevel(true);
    }

    const std::pair<float, float> point(rootBufLoc.getX(), rootBufLoc.getY());
    newSinkLocations.emplace_back(point);
    Point<double> mapKey(point.first, point.second);
    mapLocationToSink_[mapKey] = &rootBuffer;

    // Store mean LP skew target for this leaf buffer
    if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()) {
      double clusterSum = 0.0;
      for (ClockInst* ci : subInsts) {
        std::string nm = ci->getName();
        const auto sp = nm.rfind('/');
        if (sp != std::string::npos) nm = nm.substr(0, sp);
        clusterSum += cts3dDb_->getSkewTarget(nm);
      }
      clusterBufTarget_[&rootBuffer] = clusterSum / subInsts.size();
    }
    ++leafIdx;
  };

  for (const std::vector<unsigned>& cluster :
       matching.sinkClusteringSolution()) {
    if (cluster.size() == 1) {
      const std::pair<float, float>& point = points[cluster[0]];
      newSinkLocations.emplace_back(point);
      // V58 Fix A: Register singleton cluster in mapLocationToSink_.
      // Without this, the singleton's coordinate is in topLevelSinksClustered_
      // but missing from mapLocationToSink_ → computeBranchDelayTargets() and
      // getDominantTier() silently skip this sink (exact find() fails).
      // Same fix already exists in TARP path (line ~4819).
      const Point<double> sinkPt((double) point.first, (double) point.second);
      if (mapLocationToSink_.find(sinkPt) == mapLocationToSink_.end()) {
        constexpr double eps = 0.01;
        for (auto& [key, sink] : mapLocationToSink_) {
          if (std::abs(key.getX() - sinkPt.getX()) < eps
              && std::abs(key.getY() - sinkPt.getY()) < eps) {
            mapLocationToSink_[sinkPt] = sink;
            break;
          }
        }
      }
      clusterCount++;
      continue;
    }
    if (cluster.size() > 1) {
      // Collect ClockInst pointers for this cluster
      std::vector<ClockInst*> clusterClockInsts;
      for (auto point_idx : cluster) {
        const std::pair<double, double>& point = points[point_idx];
        const Point<double> mapPoint(point.first, point.second);
        if (mapLocationToSink_.find(mapPoint) == mapLocationToSink_.end()) {
          // FloatFix: epsilon fallback for float→double key mismatch
          bool found = false;
          constexpr double eps = 0.01;
          for (auto& [key, sink] : mapLocationToSink_) {
            if (std::abs(key.getX() - mapPoint.getX()) < eps
             && std::abs(key.getY() - mapPoint.getY()) < eps) {
              mapLocationToSink_[mapPoint] = sink;
              found = true;
              break;
            }
          }
          if (!found) {
            // V58: eps=1e-4 failed — find absolute nearest entry (no distance limit).
            // Skipping would leave this sink unconnected → SIGSEGV in writeDummyLoadsToDb.
            double bestDist = std::numeric_limits<double>::max();
            for (auto& [key, sink] : mapLocationToSink_) {
              double d = std::abs(key.getX() - mapPoint.getX())
                       + std::abs(key.getY() - mapPoint.getY());
              if (d < bestDist) {
                bestDist = d;
                mapLocationToSink_[mapPoint] = sink;
                found = true;
              }
            }
            if (found) {
              logger_->warn(CTS, 79,
                  "Sink at ({:.6f}, {:.6f}) — nearest match dist={:.6f}",
                  mapPoint.getX(), mapPoint.getY(), bestDist);
            }
          }
        }
        clusterClockInsts.push_back(mapLocationToSink_[mapPoint]);
      }

      // Detect cross-tier clustering
      int tier0_count = 0, tier1_count = 0;
      for (const ClockInst* inst : clusterClockInsts) {
        const int inst_tier = cts3dDb_->getInstTier(inst->getDbInst());
        if (inst_tier == 0) tier0_count++;
        else if (inst_tier == 1) tier1_count++;
      }
      if (tier0_count > 0 && tier1_count > 0) {
        logger_->warn(CTS, 571,
                      "Cross-tier cluster detected: cluster={}, tier0={}, tier1={}, "
                      "dominant_tier={}, HB_penalty={:.1f}um",
                      clusterCount, tier0_count, tier1_count,
                      cts3dDb_->getDominantTier(clusterClockInsts),
                      hb_equivalent_dist);
      }

      // Lambda for V37 target-aware split + leaf buffer
      // creation.  Applies target-spread splitting to any (sub-)cluster, then
      // creates leaf buffers.  Used by both tier-split sub-clusters and
      // unsplit clusters so that tier split no longer skips target split.
      auto targetSplitOrCreate = [&](std::vector<ClockInst*>& insts,
                                     const std::vector<unsigned>& idxs) {
        if (enableTargetSplit_ && insts.size() > 2
            && cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()) {
          struct FFTarget {
            unsigned pointIdx;
            ClockInst* inst;
            double target;
          };
          std::vector<FFTarget> ffTargets;
          ffTargets.reserve(insts.size());
          double tMin = std::numeric_limits<double>::max();
          double tMax = -std::numeric_limits<double>::max();
          for (size_t i = 0; i < insts.size(); ++i) {
            std::string nm = insts[i]->getName();
            const auto sp = nm.rfind('/');
            if (sp != std::string::npos) nm = nm.substr(0, sp);
            const double t = cts3dDb_->getSkewTarget(nm);
            ffTargets.push_back({idxs[i], insts[i], t});
            tMin = std::min(tMin, t);
            tMax = std::max(tMax, t);
          }
          const double spread = tMax - tMin;

          if (spread > splitThreshold2wayNs_) {
            std::sort(ffTargets.begin(), ffTargets.end(),
                      [](const FFTarget& a, const FFTarget& b) {
                        return a.target < b.target;
                      });

            int numSubs = 1;
            if (spread > splitThreshold3wayNs_ && ffTargets.size() >= 6) {
              numSubs = 3;
              ++splitCount_3way;
            } else if (spread > splitThreshold2wayNs_ && ffTargets.size() >= 4) {
              numSubs = 2;
              ++splitCount_2way;
            }

            if (numSubs > 1) {
              const size_t total = ffTargets.size();
              for (int s = 0; s < numSubs; ++s) {
                const size_t start = s * total / numSubs;
                const size_t end = (s + 1) * total / numSubs;
                if (start >= end) continue;

                std::vector<ClockInst*> splitInsts;
                std::vector<unsigned> splitIdxs;
                for (size_t k = start; k < end; ++k) {
                  splitInsts.push_back(ffTargets[k].inst);
                  splitIdxs.push_back(ffTargets[k].pointIdx);
                }
                if (splitInsts.size() < 2 && s + 1 < numSubs) {
                  continue;
                }

                logger_->info(CTS, 413,
                              "Sub-cluster {}.{}: {} sinks, "
                              "target=[{:.4f}, {:.4f}]ns (spread={:.4f}ns)",
                              clusterCount, s, splitInsts.size(),
                              ffTargets[start].target,
                              ffTargets[end - 1].target,
                              ffTargets[end - 1].target
                                  - ffTargets[start].target);
                createLeafBuffer(splitInsts, splitIdxs, clusterCount);
              }
              return;
            }
          }
        }
        // No target split needed: create leaf buffer directly
        createLeafBuffer(insts, idxs, clusterCount);
      };

      // Tier-aware cluster split.
      // If cluster has both bottom and upper FFs, split into two sub-clusters.
      // Each sub-cluster also gets V37 target-aware splitting if enabled.
      bool wasSplit = false;
      if (enableTierSplit_ && tier0_count > 0 && tier1_count > 0) {
        wasSplit = true;

        // Partition into bottom and upper sub-clusters
        std::vector<ClockInst*> bottomInsts, upperInsts;
        std::vector<unsigned> bottomIdxs, upperIdxs;
        for (size_t i = 0; i < cluster.size(); ++i) {
          const int inst_tier = cts3dDb_->getInstTier(
              clusterClockInsts[i]->getDbInst());
          if (inst_tier == 0) {
            bottomInsts.push_back(clusterClockInsts[i]);
            bottomIdxs.push_back(cluster[i]);
          } else {
            upperInsts.push_back(clusterClockInsts[i]);
            upperIdxs.push_back(cluster[i]);
          }
        }

        logger_->info(CTS, 377,
            "Tier split cluster {}: bottom={}, upper={}",
            clusterCount,
            bottomInsts.size(), upperInsts.size());

        // Apply V37 target split to each tier sub-cluster
        if (!bottomInsts.empty()) targetSplitOrCreate(bottomInsts, bottomIdxs);
        if (!upperInsts.empty()) targetSplitOrCreate(upperInsts, upperIdxs);
      }

      // Non-tier-split: apply V37 target split or direct leaf buffer
      if (!wasSplit) {
        targetSplitOrCreate(clusterClockInsts, cluster);
      }
    }
  }

  // Log sub-clustering statistics
  if (splitCount_2way > 0 || splitCount_3way > 0) {
    logger_->info(CTS, 414,
                  "Target-aware sub-clustering: {} 2-way splits, "
                  "{} 3-way splits, total leaf buffers = {}",
                  splitCount_2way, splitCount_3way, clusterCount);
  }

  topLevelSinksClustered_ = std::move(newSinkLocations);
  if (clusterCount) {
    treeBufLevels_++;
  }

  logger_->info(CTS,
                19,
                " Total number of sinks after clustering: {}.",
                topLevelSinksClustered_.size());

  // Verify every entry in topLevelSinksClustered_ exists in mapLocationToSink_
  int missingCount = 0;
  for (size_t si = 0; si < topLevelSinksClustered_.size(); ++si) {
    const auto& sinkPair = topLevelSinksClustered_[si];
    Point<double> key(sinkPair.first, sinkPair.second);
    if (mapLocationToSink_.find(key) == mapLocationToSink_.end()) {
      logger_->warn(CTS, 418,
                    "Debug: topLevelSinksClustered_[{}] = ({:.10f}, {:.10f}) "
                    "NOT in mapLocationToSink_ (size={})",
                    si, sinkPair.first, sinkPair.second,
                    mapLocationToSink_.size());
      ++missingCount;
    }
  }
  if (missingCount > 0) {
    logger_->warn(CTS, 419,
                  "Debug: {} / {} sinks missing from mapLocationToSink_!",
                  missingCount, topLevelSinksClustered_.size());
  } else {
    logger_->info(CTS, 427,
                  "Debug: All {} sinks verified in mapLocationToSink_ (mapSize={})",
                  topLevelSinksClustered_.size(), mapLocationToSink_.size());
  }
}

Point<double> HTreeBuilder::resolveLocationCollision(
    const Point<double>& legalCenter) const
{
  Point<double> resolvedLocation = legalCenter;

  // Collision check and jittering to ensure unique coordinates
  unsigned jitterCount = 0;
  while (true) {
    // FloatFix: removed (float) cast to match map key precision (double)
    Point<double> checkKey(resolvedLocation.getX(),
                           resolvedLocation.getY());
    if (mapLocationToSink_.find(checkKey) == mapLocationToSink_.end()) {
      break;
    }
    jitterCount++;
    double offset = (double) jitterCount / wireSegmentUnit_;
    resolvedLocation.setX(legalCenter.getX() + offset);
  }

  if (jitterCount > 0) {
    debugPrint(
        logger_,
        utl::CTS,
        "clustering",
        1,
        "Resolved collision via jittering ({} times): ({}, {}) -> ({}, {})",
        jitterCount,
        legalCenter.getX(),
        legalCenter.getY(),
        resolvedLocation.getX(),
        resolvedLocation.getY());
  }
  return resolvedLocation;
}

// Removed getDominantTierFromInsts, getDominantTierFromSinkLocs,
// getDominantTierFromClockSinks — moved to Cts3DDatabase::getDominantTier() and
// Cts3DDatabase::getDominantTierFromClock()

void HTreeBuilder::initSinkRegion()
{
  const unsigned wireSegmentUnitInDbu = techChar_->getLengthUnit();
  const int dbUnits = options_->getDbUnits();
  wireSegmentUnit_ = wireSegmentUnitInDbu;

  double clusterDiameter = (type_ == TreeType::MacroTree)
                               ? options_->getMacroMaxDiameter()
                               : options_->getMaxDiameter();
  unsigned clusterSize = (type_ == TreeType::MacroTree)
                             ? options_->getMacroSinkClusteringSize()
                             : options_->getSinkClusteringSize();
  unsigned min_clustering_sinks = (type_ == TreeType::MacroTree)
                                      ? min_clustering_macro_sinks_
                                      : min_clustering_sinks_;

  logger_->info(CTS,
                20,
                " Wire segment unit: {}  dbu ({} um).",
                wireSegmentUnit_,
                wireSegmentUnitInDbu / dbUnits);

  if (options_->isSimpleSegmentEnabled()) {
    const int remainingLength
        = options_->getBufferDistance() / (wireSegmentUnitInDbu * 2);
    logger_->info(CTS,
                  21,
                  " Distance between buffers: {} units ({} um).",
                  remainingLength,
                  static_cast<int>(options_->getBufferDistance() / dbUnits));
    if (options_->isVertexBuffersEnabled()) {
      const int vertexBufferLength
          = options_->getVertexBufferDistance() / (wireSegmentUnitInDbu * 2);
      logger_->info(
          CTS,
          22,
          " Branch length for Vertex Buffer: {} units ({} um).",
          vertexBufferLength,
          static_cast<int>(options_->getVertexBufferDistance() / dbUnits));
    }
  }

  std::vector<std::pair<float, float>> topLevelSinks;
  std::vector<const ClockInst*> sinkInsts;
  initTopLevelSinks(topLevelSinks, sinkInsts);

  const float maxDiameter = (clusterDiameter * dbUnits) / wireSegmentUnit_;
  // clang-format off
  debugPrint(logger_, CTS, "clustering", 1, "maxDiameter={:0.3f} = "
             "origMaxDiam={} * dbUnits={} / wireSegmentUnit_={}",
             maxDiameter, clusterDiameter, dbUnits,
             wireSegmentUnit_);
  // clang-format on

  preSinkClustering(topLevelSinks, sinkInsts, maxDiameter, clusterSize);
  if (topLevelSinks.size() <= min_clustering_sinks
      || !(options_->getSinkClustering())) {
    Box<int> sinkRegionDbu = clock_.computeSinkRegion();
    logger_->info(CTS, 23, " Original sink region: {}.", sinkRegionDbu);

    sinkRegion_ = sinkRegionDbu.normalize(1.0 / wireSegmentUnit_);
  } else {
    if (topLevelSinksClustered_.size() > 400
        && options_->getSinkClusteringLevels() > 0) {
      std::vector<std::pair<float, float>> secondLevelLocs;
      std::vector<const ClockInst*> secondLevelInsts;
      initSecondLevelSinks(secondLevelLocs, secondLevelInsts);
      preSinkClustering(secondLevelLocs,
                        secondLevelInsts,
                        maxDiameter * 4,
                        std::ceil(std::sqrt(clusterSize)),
                        true);
    }
    sinkRegion_ = clock_.computeSinkRegionClustered(topLevelSinksClustered_);
  }
  logger_->info(CTS, 24, " Normalized sink region: {}.", sinkRegion_);
  logger_->info(CTS, 25, "    Width:  {:.4f}.", sinkRegion_.getWidth());
  logger_->info(CTS, 26, "    Height: {:.4f}.", sinkRegion_.getHeight());
}

static void plotBlockage(std::ofstream& file,
                         odb::dbDatabase* db_,
                         int scalingFactor)
{
  unsigned i = 0;
  for (odb::dbBlockage* blockage : db_->getChip()->getBlock()->getBlockages()) {
    odb::dbBox* bbox = blockage->getBBox();
    int x = bbox->xMin() / scalingFactor;
    int y = bbox->yMin() / scalingFactor;
    int w = bbox->xMax() / scalingFactor - bbox->xMin() / scalingFactor;
    int h = bbox->yMax() / scalingFactor - bbox->yMin() / scalingFactor;
    file << i++ << " " << x << " " << y << " " << w << " " << h
         << " block  scalingFactor=";
    file << scalingFactor << " " << blockage->getId() << '\n';
  }
}

// distance to move sinks from old loc to new loc
double HTreeBuilder::weightedDistance(const Point<double>& newLoc,
                                      const Point<double>& oldLoc,
                                      const std::vector<Point<double>>& sinks)
{
  double dist = 0;
  for (const Point<double>& sink : sinks) {
    dist += computeDist(newLoc, sink);
    dist += computeDist(newLoc, oldLoc);
  }
  return dist;
}

static void plotSinks(std::ofstream& file,
                      const std::vector<Point<double>>& sinks)
{
  unsigned cnt = 0;
  for (const Point<double>& pt : sinks) {
    double x = pt.getX();
    double y = pt.getY();
    double w = 1;
    double h = 1;
    auto name = "sink_";
    file << cnt++ << " " << x << " " << y << " " << w << " " << h;
    file << " " << name << '\n';
  }
}

unsigned HTreeBuilder::findSibling(LevelTopology& topology,
                                   unsigned i,
                                   unsigned par)
{
  for (unsigned idx = 0; idx < topology.getBranchingPointSize(); ++idx) {
    unsigned k = topology.getBranchingPointParentIdx(idx);
    if (idx != i && k == par) {
      return idx;
    }
  }
  return i;
}

void HTreeBuilder::scalePosition(Point<double>& loc,
                                 const Point<double>& parLoc,
                                 double leng,
                                 double scale)
{
  double px = parLoc.getX();
  double py = parLoc.getY();
  double ax = loc.getX();
  double ay = loc.getY();

  double d = computeDist(loc, parLoc);
  double x, y;
  if (d > 0) {  // yy8
    double delta = d * scale;
    double dx = ax - px;
    double dy = ay - py;
    dx += (dx > 0) ? delta : -delta;
    dy += (dy > 0) ? -delta : delta;
    double scale = leng / d;
    x = px + dx * scale;
    y = py + dy * scale;
  } else {
    x = px + leng / 2;
    y = py + leng / 2;
  }
  loc.setX(x);
  loc.setY(y);
}

static void setSiblingPosition(const Point<double>& a,
                               Point<double>& b,
                               const Point<double>& parLoc)
{
  double px = parLoc.getX();
  double py = parLoc.getY();
  double ax = a.getX();
  double ay = a.getY();
  double bx = 2 * px - ax;
  double by = 2 * py - ay;
  b.setX(bx);
  b.setY(by);
}

// Balance the two branches on the very top level
void HTreeBuilder::adjustToplevelTopology(Point<double>& a,
                                          Point<double>& b,
                                          const Point<double>& parLoc)
{
  double da = computeDist(a, parLoc);
  double db = computeDist(b, parLoc);
  if (da < db) {
    setSiblingPosition(a, b, parLoc);
  } else {
    setSiblingPosition(b, a, parLoc);
  }
}

void HTreeBuilder::findLegalLocations(const Point<double>& parentPoint,
                                      const Point<double>& branchPoint,
                                      double x1,
                                      double y1,
                                      double x2,
                                      double y2,
                                      std::vector<Point<double>>& points)
{
  // add 4 corners of blockage
  addCandidateLoc(x1, y1, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc(x1, y2, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc(x2, y2, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc(x2, y1, parentPoint, x1, y1, x2, y2, points);

  // add straightline neighbors
  double bx = branchPoint.getX();
  double by = branchPoint.getY();
  addCandidateLoc(bx, y1, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc(bx, y2, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc(x1, by, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc(x2, by, parentPoint, x1, y1, x2, y2, points);

  double px = parentPoint.getX();
  double py = parentPoint.getY();
  double dx = px - bx;
  double dy = py - by;
  double m = dy / dx;
  // y = m*(x-bx) + by
  addCandidateLoc(x1, m * (x1 - bx) + by, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc(x2, m * (x2 - bx) + by, parentPoint, x1, y1, x2, y2, points);
  // x = (y-by)/m + bx
  addCandidateLoc((y1 - by) / m + bx, y1, parentPoint, x1, y1, x2, y2, points);
  addCandidateLoc((y2 - by) / m + bx, y2, parentPoint, x1, y1, x2, y2, points);
  // clang-format off
  if (logger_->debugCheck(utl::CTS, "legalizer", 3)) {
    logger_->report("    branchPt:{} is not legal, parentPt:{} blockages:({:0.3f} {:0.3f}) "
        "({:0.3f} {:0.3f})", branchPoint, parentPoint, x1, y1, x2, y2);
    for (Point<double> point : points) {
      logger_->report("      FLL candiate {}", point);
    }
    // clang-format on
  }
}

Point<double> HTreeBuilder::findBestLegalLocation(
    double targetDist,
    const Point<double>& branchPoint,
    const Point<double>& parentPoint,
    const std::vector<Point<double>>& legalLocations,
    const std::vector<Point<double>>& sinks,
    double x1,
    double y1,
    double x2,
    double y2,
    int scalingFactor,
    odb::Direction2D direction)
{
  Point<double> best(0.0, 0.0);
  double minDiff = std::numeric_limits<double>::max();
  for (const Point<double>& loc : legalLocations) {
    double dist = computeDist(loc, parentPoint);
    dist += computeDist(loc, branchPoint);
    double diff = abs(dist - targetDist);
    // clang-format off
    if (logger_->debugCheck(utl::CTS, "legalizer", 3)) {
      logger_->report("      Loc {}: curr dist={:0.3f} target dist={:0.3f} sink"
		      " dist={:0.3f}", loc, dist, targetDist,
		      weightedDistance(loc, branchPoint, sinks));
    }
    // clang-format on
    if (diff < minDiff) {
      minDiff = diff;
      best = loc;
    }
  }

  return adjustBestLegalLocation(targetDist,
                                 best,
                                 parentPoint,
                                 sinks,
                                 x1,
                                 y1,
                                 x2,
                                 y2,
                                 scalingFactor,
                                 direction);
}

// Adjust buffer location in two steps:
// step1: try moving buffer to existing blockage boundary (less expensive)
// step2: try moving buffer beyond existing blockage boundary (more expensive)
// In both steps, the first priority is to match target distance from current
// point to parent point.  The second priority is to lower the weighted sink
// distance
Point<double> HTreeBuilder::adjustBestLegalLocation(
    double targetDist,
    const Point<double>& currLoc,
    const Point<double>& parentPoint,
    const std::vector<Point<double>>& sinks,
    double x1,
    double y1,
    double x2,
    double y2,
    int scalingFactor,
    odb::Direction2D direction)
{
  if (fuzzyEqual(targetDist, computeDist(currLoc, parentPoint))) {
    return currLoc;
  }

  // try moving along blockage boundary
  Point<double> bestLoc = currLoc;
  if (adjustAlongBlockage(targetDist,
                          currLoc,
                          parentPoint,
                          sinks,
                          x1,
                          y1,
                          x2,
                          y2,
                          scalingFactor,
                          bestLoc)) {
    return bestLoc;
  }

  // try moving beyond blockage boundary
  return adjustBeyondBlockage(
      currLoc, parentPoint, targetDist, sinks, scalingFactor, direction);
}

bool HTreeBuilder::adjustAlongBlockage(double targetDist,
                                       const Point<double>& currLoc,
                                       const Point<double>& parentPoint,
                                       const std::vector<Point<double>>& sinks,
                                       double x1,
                                       double y1,
                                       double x2,
                                       double y2,
                                       int scalingFactor,
                                       Point<double>& bestLoc)
{
  Point<double> newLoc = currLoc;
  // clang-format off
  debugPrint(logger_, CTS, "legalizer", 3, "{} currDist={:0.3f} != "
	     "targetDist={:0.3f}, adjustAlongBlockage...", currLoc,
	     computeDist(currLoc, parentPoint), targetDist);
  // clang-format on
  double x = currLoc.getX();
  double y = currLoc.getY();
  double px = parentPoint.getX();
  double py = parentPoint.getY();

  Point<double> noBest(std::numeric_limits<double>::min(),
                       std::numeric_limits<double>::min());
  bestLoc = noBest;
  double bestSinkDist = std::numeric_limits<double>::max();
  double sinkDist = 0.0;
  std::vector<Point<double>> candidates;

  // move along y axis within blockage
  double newX = x;
  double newY = py - targetDist + abs(px - newX);
  newLoc.setX(newX);
  newLoc.setY(newY);
  candidates.emplace_back(newLoc);  // trial y#1

  newY = py + targetDist - abs(px - newX);
  newLoc.setY(newY);
  candidates.emplace_back(newLoc);  // trial y#2

  // move along x axis within blockage
  newY = y;
  newX = px - targetDist + abs(py - newY);
  newLoc.setX(newX);
  newLoc.setY(newY);
  candidates.emplace_back(newLoc);  // trial x#1

  newX = px + targetDist - abs(py - newY);
  newLoc.setX(newX);
  candidates.emplace_back(newLoc);  // trial x#2

  for (const Point<double>& candidate : candidates) {
    checkLegalityAndCostSpecial(currLoc,
                                candidate,
                                parentPoint,
                                targetDist,
                                sinks,
                                scalingFactor,
                                x1,
                                y1,
                                x2,
                                y2,
                                bestLoc,
                                sinkDist,
                                bestSinkDist);
  }

  return (bestLoc != noBest);
}

void HTreeBuilder::checkLegalityAndCostSpecial(
    const Point<double>& oldLoc,
    const Point<double>& newLoc,
    const Point<double>& parentPoint,
    double targetDist,
    const std::vector<Point<double>>& sinks,
    int scalingFactor,
    double x1,
    double y1,
    double x2,
    double y2,
    Point<double>& bestLoc,
    double& sinkDist,
    double& bestSinkDist)
{
  if (fuzzyEqual(computeDist(newLoc, parentPoint), targetDist)
      && checkLegalitySpecial(newLoc, x1, y1, x2, y2, scalingFactor)) {
    sinkDist = weightedDistance(newLoc, oldLoc, sinks);
    if (sinkDist < bestSinkDist) {
      bestLoc = newLoc;
      bestSinkDist = sinkDist;
    }
    // clang-format off
    debugPrint(logger_, CTS, "legalizer", 3, "adjustBestLegalLoc: branchPt "
	       "move:{}=>{} is legal, dist={:0.3f}, sinkDist={:0.3f}",
	       oldLoc, newLoc, targetDist, sinkDist);
    // clang-format on
  }
  // clang-format off
  debugPrint(logger_, CTS, "legalizer", 3, "adjustBestLegalLoc: branchPt "
	     "move:{}=>{} is illegal or dist {:0.3f} != {:0.3f}",
	     oldLoc, newLoc, computeDist(newLoc, parentPoint), targetDist);
  // clang-format on
}

// 1) Branch point couldn't be legalized by simply moving it along blockage
//    boundary, or
// 2) Branch point is legal but parent point has moved, so it is necessary to
//    move it to match target topology length
// In both cases, move branch point to match target distance and minimize
// weighted sink distance
Point<double> HTreeBuilder::adjustBeyondBlockage(
    const Point<double>& branchPoint,
    const Point<double>& parentPoint,
    double targetDist,
    const std::vector<Point<double>>& sinks,
    int scalingFactor,
    odb::Direction2D direction)
{
  double px = parentPoint.getX();
  double py = parentPoint.getY();
  std::vector<Point<double>> candidates;
  Point<double> point = branchPoint;

  // try points that are on edges of "Manhattan square"
  // with the parent point in the center
  //              p16 p1 p9
  //            p8         p5
  //          p15            p10
  //        p4        pp       p2
  //         p14             p11
  //           p7           p6
  //              p13 p3 p12
  double leng50 = targetDist * 0.5;
  double leng25 = targetDist * 0.25;
  double leng75 = targetDist * 0.75;
  switch (direction) {
    case odb::Direction2D::North:
      addCandidatePoint(px, py + targetDist, point, candidates);       // p1
      addCandidatePoint(px + leng50, py + leng50, point, candidates);  // p5
      addCandidatePoint(px - leng50, py + leng50, point, candidates);  // p8
      addCandidatePoint(px + leng25, py + leng75, point, candidates);  // p9
      addCandidatePoint(px - leng25, py + leng75, point, candidates);  // p16
      break;
    case odb::Direction2D::East:
      addCandidatePoint(px + targetDist, py, point, candidates);       // p2
      addCandidatePoint(px + leng50, py + leng50, point, candidates);  // p5
      addCandidatePoint(px + leng50, py - leng50, point, candidates);  // p6
      addCandidatePoint(px + leng75, py + leng25, point, candidates);  // p10
      addCandidatePoint(px + leng75, py - leng25, point, candidates);  // p11
      break;
    case odb::Direction2D::South:
      addCandidatePoint(px, py - targetDist, point, candidates);       // p3
      addCandidatePoint(px + leng50, py - leng50, point, candidates);  // p6
      addCandidatePoint(px - leng50, py - leng50, point, candidates);  // p7
      addCandidatePoint(px + leng25, py - leng75, point, candidates);  // p12
      addCandidatePoint(px - leng25, py - leng75, point, candidates);  // p13
      break;
    default:
      addCandidatePoint(px - targetDist, py, point, candidates);       // p4
      addCandidatePoint(px - leng50, py - leng50, point, candidates);  // p7
      addCandidatePoint(px - leng50, py + leng50, point, candidates);  // p8
      addCandidatePoint(px - leng75, py - leng25, point, candidates);  // p14
      addCandidatePoint(px - leng75, py + leng25, point, candidates);  // p15
      // odb::Direction2D::West
      break;
  }

  // check if any corners of Manhanttan square are inside some blockage
  // if so add candidate points that intersect blockage and Manhattan square
  // p1 is top corner

  if (direction == odb::Direction2D(odb::Direction2D::North)) {
    point.setX(px);
    point.setY(py + targetDist);
    addCandidatePointsAlongBlockage(point,
                                    parentPoint,
                                    targetDist,
                                    scalingFactor,
                                    candidates,
                                    odb::Direction2D::North);
  }

  // p2 is right corner
  if (direction == odb::Direction2D(odb::Direction2D::East)) {
    point.setX(px + targetDist);
    point.setY(py);
    addCandidatePointsAlongBlockage(point,
                                    parentPoint,
                                    targetDist,
                                    scalingFactor,
                                    candidates,
                                    odb::Direction2D::East);
  }

  // p3 is bottom corner
  if (direction == odb::Direction2D(odb::Direction2D::South)) {
    point.setX(px);
    point.setY(py - targetDist);
    addCandidatePointsAlongBlockage(point,
                                    parentPoint,
                                    targetDist,
                                    scalingFactor,
                                    candidates,
                                    odb::Direction2D::South);
  }

  // p4 is left corner
  if (direction == odb::Direction2D(odb::Direction2D::West)) {
    point.setX(px - targetDist);
    point.setY(py);
    addCandidatePointsAlongBlockage(point,
                                    parentPoint,
                                    targetDist,
                                    scalingFactor,
                                    candidates,
                                    odb::Direction2D::West);
  }

  // try moving cell along x or y, with some offset
  double bx = branchPoint.getX();
  double by = branchPoint.getY();
  Point<double> noBest(std::numeric_limits<double>::min(),
                       std::numeric_limits<double>::min());
  Point<double> bestLoc = noBest;
  double bestSinkDist = std::numeric_limits<double>::max();
  double sinkDist = 0.0;
  // get small x and y offset to find "channel" through blockages
  double minX = 10.0 * getBufferWidth();
  double minY = 10.0 * getBufferHeight();

  // try small offset in x direction
  Point<double> newLoc(bx, by);
  double newX = bx + minX;
  double newY = py - targetDist + abs(px - newX);
  newLoc.setX(newX);
  newLoc.setY(newY);
  if (direction == odb::Direction2D(odb::Direction2D::South)) {
    candidates.emplace_back(newLoc);  // trial x#1
  }

  newLoc.setY(py + targetDist - abs(px - newX));
  if (direction == odb::Direction2D(odb::Direction2D::North)) {
    candidates.emplace_back(newLoc);  // trial x#2
  }

  newX = bx - minX;
  newY = py - targetDist + abs(px - newX);
  newLoc.setX(newX);
  newLoc.setY(newY);
  if (direction == odb::Direction2D(odb::Direction2D::South)) {
    candidates.emplace_back(newLoc);  // trial x#3
  }

  newLoc.setY(py + targetDist - abs(px - newX));
  if (direction == odb::Direction2D(odb::Direction2D::North)) {
    candidates.emplace_back(newLoc);  // trial x#4
  }

  // try small offset in y direction
  newY = by + minY;
  newX = px - targetDist + abs(py - newY);
  newLoc.setX(newX);
  newLoc.setY(newY);
  if (direction == odb::Direction2D(odb::Direction2D::West)) {
    candidates.emplace_back(newLoc);  // trial y#1
  }

  newLoc.setX(px + targetDist - abs(py - newY));
  if (direction == odb::Direction2D(odb::Direction2D::East)) {
    candidates.emplace_back(newLoc);  // trial y#2
  }

  newY = by - minY;
  newX = px - targetDist + abs(py - newY);
  newLoc.setX(newX);
  newLoc.setY(newY);
  if (direction == odb::Direction2D(odb::Direction2D::West)) {
    candidates.emplace_back(newLoc);  // trial y#3
  }

  newLoc.setX(px + targetDist - abs(py - newY));
  if (direction == odb::Direction2D(odb::Direction2D::East)) {
    candidates.emplace_back(newLoc);  // trial y#4
  }

  for (const Point<double>& candidate : candidates) {
    checkLegalityAndCost(branchPoint,
                         candidate,
                         parentPoint,
                         targetDist,
                         sinks,
                         scalingFactor,
                         bestLoc,
                         sinkDist,
                         bestSinkDist);
  }

  if (bestLoc != noBest) {
    return bestLoc;
  }
  return branchPoint;
}

void HTreeBuilder::addCandidatePointsAlongBlockage(
    const Point<double>& point,
    const Point<double>& parentPoint,
    double targetDist,
    int scalingFactor,
    std::vector<Point<double>>& candidates,
    odb::Direction2D direction)
{
  double x1, y1, x2, y2;
  if (findBlockage(point, scalingFactor, x1, y1, x2, y2)) {
    Point<double> point2 = point;
    double px = parentPoint.getX();
    double py = parentPoint.getY();
    // clang-format off
    debugPrint(logger_, CTS, "legalizer", 3, "  {} corner {} of Manhattan "
	       "sqaure is inside blockage ({:0.3f} {:0.3f}) ({:0.3f} {:0.3f})",
	       direction, point, x1, y1, x2, y2);
    // clang-format on
    switch (direction) {
      case odb::Direction2D::North:
        //        ----------
        //        |        |
        //        |   top  |
        //        |   / \  |
        // (x1, y1)--x---x-(x2, y1)
        //          /  .
        //          (px, py)
        //
        // (px - x) + (y1 - py) = targetDist
        // (x - px) + (y1 - py) = targetDist
        point2.setX(-targetDist - py + px + y1);
        point2.setY(y1);
        candidates.emplace_back(point2);
        point2.setY(targetDist + py + px - y1);
        candidates.emplace_back(point2);
        break;
      case odb::Direction2D::East:
        // (x1, y2)---------
        //       \ |        |
        //         x        |
        //         |\       |
        //   .     | right  |
        // (px,py)| /       |
        //         x        |
        // (x1, y1)---------
        // (x1 - px) + (y - py) = targetDist
        // (x1 - px) + (py - y) = targetDist
        point2.setX(x1);
        point2.setY(targetDist + py + px - x1);
        candidates.emplace_back(point2);
        point2.setY(-targetDist + py - px + x1);
        candidates.emplace_back(point2);
        break;
      case odb::Direction2D::South:
        point2.setX(-targetDist + py + px - y2);
        point2.setY(y2);
        candidates.emplace_back(point2);
        point2.setX(targetDist - py + px + y2);
        candidates.emplace_back(point2);
        break;
      default:
        // odb::Direction2D::West
        point2.setX(x2);
        point2.setY(targetDist + py - px + x2);
        candidates.emplace_back(point2);
        point2.setY(-targetDist + py + px - x2);
        candidates.emplace_back(point2);
        break;
    }
  }
}

void HTreeBuilder::checkLegalityAndCost(const Point<double>& oldLoc,
                                        const Point<double>& newLoc,
                                        const Point<double>& parentPoint,
                                        double targetDist,
                                        const std::vector<Point<double>>& sinks,
                                        int scalingFactor,
                                        Point<double>& bestLoc,
                                        double& sinkDist,
                                        double& bestSinkDist)
{
  if (fuzzyEqual(computeDist(newLoc, parentPoint), targetDist)
      && checkLegalityLoc(newLoc, scalingFactor)) {
    sinkDist = weightedDistance(newLoc, oldLoc, sinks);
    if (sinkDist < bestSinkDist) {
      bestLoc = newLoc;
      bestSinkDist = sinkDist;
    }
    // clang-format off
    debugPrint(logger_, CTS, "legalizer", 3, "adjustBeyondBlockage: branchPt "
	       "move:{}=>{} is legal, dist={:0.3f}, sinkDist={:0.3f}",
	       oldLoc, newLoc, targetDist, sinkDist);
  } else {
    debugPrint(logger_, CTS, "legalizer", 3, "adjustBeyondBlockage: branchPt "
	       "move:{}=>{} is illegal or dist {:0.3f} != {:0.3f}",
	       oldLoc, newLoc, computeDist(newLoc, parentPoint), targetDist);
    // clang-format on
  }
}

void HTreeBuilder::legalizeDummy()
{
  Point<double> topLevelBufferLoc = sinkRegion_.getCenter();
  for (int levelIdx = 0; levelIdx < topologyForEachLevel_.size(); ++levelIdx) {
    LevelTopology& topology = topologyForEachLevel_[levelIdx];

    for (unsigned idx = 0; idx < topology.getBranchingPointSize(); ++idx) {
      Point<double>& branchPoint = topology.getBranchingPoint(idx);
      unsigned parentIdx = topology.getBranchingPointParentIdx(idx);

      // clang-format off
      Point<double> parentPoint
          = (levelIdx == 0)
                ? topLevelBufferLoc
                : topologyForEachLevel_[levelIdx - 1].getBranchingPoint(
                      parentIdx);
      // clang-format on

      const std::vector<Point<double>>& sinks
          = topology.getBranchSinksLocations(idx);

      double leng = topology.getLength();
      Point<double>& sibLoc = findSiblingLoc(topology, idx, parentIdx);

      double d1 = computeDist(branchPoint, sibLoc);
      double d2 = computeDist(branchPoint, parentPoint);
      bool overlap = d1 == 0 || d2 == 0;
      bool dummy = sinks.empty();  // dummy buffers drive no sinks

      // not important, can be removed later?
      if (dummy) {
        setSiblingPosition(sibLoc, branchPoint, parentPoint);
        scalePosition(branchPoint, parentPoint, leng, 0.1);
      } else if (overlap) {
        scalePosition(branchPoint, parentPoint, leng, 0.1);
      } else {
        continue;
      }

      double x1, y1, x2, y2;
      int scalingFactor = wireSegmentUnit_;
      if (!isOccupiedLoc(branchPoint)
          && findBlockage(branchPoint, scalingFactor, x1, y1, x2, y2)) {
        Point<double> legalBranchPoint(branchPoint);
        std::vector<Point<double>> legalLocations;
        findLegalLocations(
            parentPoint, branchPoint, x1, y1, x2, y2, legalLocations);
        legalBranchPoint = findBestLegalLocation(topology.getLength(),
                                                 branchPoint,
                                                 parentPoint,
                                                 legalLocations,
                                                 sinks,
                                                 x1,
                                                 y2,
                                                 x2,
                                                 y2,
                                                 scalingFactor,
                                                 odb::Direction2D::North);
        double d = computeDist(legalBranchPoint, parentPoint);
        // clang-format off
        debugPrint(logger_, CTS, "legalizer", 1,
            "legalizeDummy level index {}: {}->{} d={:0.3f}, leng={:0.3f},"
		   "ratio={:0.3f}", levelIdx, branchPoint, legalBranchPoint,
		   d, leng, d / leng);
        // clang-format on
        commitMoveLoc(branchPoint, legalBranchPoint);
        branchPoint.setX(legalBranchPoint.getX());
        branchPoint.setY(legalBranchPoint.getY());
      }
    }
  }
}

void HTreeBuilder::legalize()
{
  if (logger_->debugCheck(utl::CTS, "legalizer", 3)) {
    logger_->report("HTree before legalization -------");
    printHTree();
  }
  // make sure top level buffer is legal
  Point<double> oldTopBufferLoc = sinkRegion_.getCenter();
  Point<double> newTopBufferLoc
      = legalizeOneBuffer(oldTopBufferLoc, options_->getRootBuffer());
  sinkRegion_.setCenter(newTopBufferLoc);
  commitMoveLoc(oldTopBufferLoc, newTopBufferLoc);
  // clang-format off
  debugPrint(logger_, CTS, "legalizer", 3, "legalize: top buf loc:{}->{}",
	     oldTopBufferLoc, newTopBufferLoc);
  // clang-format on
  for (int levelIdx = 0; levelIdx < topologyForEachLevel_.size(); ++levelIdx) {
    LevelTopology& topology = topologyForEachLevel_[levelIdx];

    for (unsigned bufferIdx = 0; bufferIdx < topology.getBranchingPointSize();
         ++bufferIdx) {
      // bufferIdx is the buffer id at level levelIdx
      Point<double>& branchPoint = topology.getBranchingPoint(bufferIdx);
      unsigned parentIdx = topology.getBranchingPointParentIdx(bufferIdx);

      // clang-format off
      Point<double> parentPoint
          = (levelIdx == 0)
                ? newTopBufferLoc
                : topologyForEachLevel_[levelIdx - 1].getBranchingPoint(
                      parentIdx);
      // clang-format on

      odb::Direction2D::Value branch_point_dir;
      if (isHorizontal(levelIdx + 1)) {
        if (branchPoint.getX() - parentPoint.getX() > 0) {
          branch_point_dir = odb::Direction2D::East;
        } else {
          branch_point_dir = odb::Direction2D::West;
        }
      } else {
        if (branchPoint.getY() - parentPoint.getY() > 0) {
          branch_point_dir = odb::Direction2D::North;
        } else {
          branch_point_dir = odb::Direction2D::South;
        }
      }

      const std::vector<Point<double>>& sinks
          = topology.getBranchSinksLocations(bufferIdx);

      double leng = computeDist(branchPoint, parentPoint);
      // clang-format off
      if (logger_->debugCheck(utl::CTS, "legalizer", 3)) {
        logger_->report("  HTree level*{}* bufId*{}*, parent:{}, branch:{}, "
			"leng:{:0.3f}, sinks:{}", levelIdx, bufferIdx,
			parentPoint, branchPoint, leng, sinks.size());
      }
      // clang-format on
      int scalingFactor = wireSegmentUnit_;
      double x1, y1, x2, y2;
      if (!isOccupiedLoc(branchPoint)
          && findBlockage(branchPoint, scalingFactor, x1, y1, x2, y2)) {
        Point<double> legalBranchPoint(branchPoint);
        std::vector<Point<double>> legalLocations;
        // find all the possible locations off the blockage
        findLegalLocations(
            parentPoint, branchPoint, x1, y1, x2, y2, legalLocations);
        // choose the best new location based on desired topology length
        legalBranchPoint = findBestLegalLocation(topology.getLength(),
                                                 branchPoint,
                                                 parentPoint,
                                                 legalLocations,
                                                 sinks,
                                                 x1,
                                                 y1,
                                                 x2,
                                                 y2,
                                                 scalingFactor,
                                                 branch_point_dir);
        // clang-format off
	debugPrint(logger_, CTS, "legalizer", 1,
		   "findBestLegalLocation branchPt:{}=>{} parentPt:{} new "
		   "branchPt is {} blockage", branchPoint, legalBranchPoint,
		   parentPoint,
		   isInsideBbox(legalBranchPoint.getX(), legalBranchPoint.getY(),
				x1, y1, x2, y2)? "inside" : "outside");
        // clang-format on
        // update branchPoint
        commitMoveLoc(branchPoint, legalBranchPoint);
        branchPoint = legalBranchPoint;
      } else if (isOccupiedLoc(branchPoint)
                 || !fuzzyEqual(leng, topology.getLength(), 0.01)) {
        // legal branch point needs adjustment if parent point moved in previous
        // level
        Point<double> newLocation(branchPoint);
        newLocation = adjustBeyondBlockage(branchPoint,
                                           parentPoint,
                                           topology.getLength(),
                                           sinks,
                                           scalingFactor,
                                           branch_point_dir);
        // clang-format off
	debugPrint(logger_, CTS, "legalizer", 3,
		   "adjustBeyondBlockage applied to legal branchPt:"
		   "{}=>{} parentPt:{} newDist={:0.3f}", branchPoint, newLocation,
		   parentPoint, computeDist(newLocation, parentPoint));
        // clang-format on
        commitMoveLoc(branchPoint, newLocation);
        branchPoint = newLocation;
      } else {
        commitLoc(branchPoint);
      }
    }
  }

  // "further" optimize the location of the "dummy" buffers that drive
  // no sinks (still needed?)
  legalizeDummy();

  if (logger_->debugCheck(utl::CTS, "legalizer", 3)) {
    logger_->report("HTree after legalization -------");
    printHTree();
  }
}

void HTreeBuilder::run()
{
  double clusterDiameter = (type_ == TreeType::MacroTree)
                               ? options_->getMacroMaxDiameter()
                               : options_->getMaxDiameter();
  unsigned clusterSize = (type_ == TreeType::MacroTree)
                             ? options_->getMacroSinkClusteringSize()
                             : options_->getSinkClusteringSize();
  bool useMaxCap = (type_ == TreeType::MacroTree)
                       ? false
                       : options_->getSinkClusteringUseMaxCap();

  logger_->info(
      CTS, 27, "Generating H-Tree topology for net {}.", clock_.getName());
  logger_->info(CTS, 28, " Total number of sinks: {}.", clock_.getNumSinks());
  if (options_->getSinkClustering()) {
    if (useMaxCap) {
      logger_->info(
          CTS, 90, " Sinks will be clustered based on buffer max cap.");
    } else {
      logger_->info(
          CTS,
          29,
          " {} sinks will be clustered in groups of up to {} and with "
          "maximum cluster diameter of {:.1f} um.",
          type_ == TreeType::MacroTree ? "Macro " : "Register",
          clusterSize,
          clusterDiameter);
    }
  }
  logger_->info(
      CTS, 30, " Number of static layers: {}.", options_->getNumStaticLayers());

  clockTreeMaxDepth_ = options_->getClockTreeMaxDepth();
  minInputCap_ = techChar_->getActualMinInputCap();
  numMaxLeafSinks_ = options_->getNumMaxLeafSinks();
  minLengthSinkRegion_ = techChar_->getMinSegmentLength() * 2;

  // skewTargetBeta_ / CTS_SKEW_TARGET_BETA removed.
  // Was dead code — computeDist() called with candidate branching points
  // not in mapLocationToSink_ → beta penalty never fired.
  // Target-aware clustering: CTS_CKMEANS_TARGET_BETA (read at line ~1725).
  if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()) {
    logger_->info(CTS, 400,
                  "{} skew targets loaded",
                  cts3dDb_->getSkewTargetCount());
  }

  // Initialize useful-skew wire adjustment parameters
  // CTS_WIRE_SKEW_SCALE: 0 = disabled (zero-skew), 1.0 = full LP-guided wire shift
  // CTS_WIRE_DELAY_PS_UM: wire delay per unit length (ps/um)
  //   Used to convert LP target difference (ns) → wire length shift (um)
  //   Reference: Fishburn 1990, LP-SAFETY clock skew optimization
  wireSkewScale_    = 0.0;
  wireDelayPerUnit_ = 0.001;  // default: 1.0 ps/um in ns/um
  if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()) {
    const char* wss_env  = std::getenv("CTS_WIRE_SKEW_SCALE");
    const char* wdpu_env = std::getenv("CTS_WIRE_DELAY_PS_UM");
    if (wss_env != nullptr) {
      wireSkewScale_ = std::atof(wss_env);
    } else {
      wireSkewScale_ = 0.0;  // V58 Fix D: default OFF — hold regression (V51b confirmed)
    }
    if (wdpu_env != nullptr) {
      wireDelayPerUnit_ = std::atof(wdpu_env) * 0.001;  // ps/um → ns/um
    }
    if (wireSkewScale_ > 0.0) {
      logger_->info(CTS, 405,
                    "Useful-skew wire adjustment enabled: scale={:.2f}, "
                    "wireDelay={:.3f}ps/um",
                    wireSkewScale_, wireDelayPerUnit_ * 1000.0);
    }
  }

  // Read target-split env vars BEFORE initSinkRegion(),
  // because preSinkClustering() (called inside initSinkRegion()) uses
  // enableTargetSplit_ at line 328.  Previously read at line 1576 — too late.
  if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()) {
    if (const char* e = std::getenv("CTS_ENABLE_TARGET_SPLIT"))
      enableTargetSplit_ = (std::atoi(e) != 0);
    if (const char* e = std::getenv("CTS_SPLIT_THRESHOLD_2WAY_PS"))
      splitThreshold2wayNs_ = std::atof(e) * 0.001;  // ps → ns
    if (const char* e = std::getenv("CTS_SPLIT_THRESHOLD_3WAY_PS"))
      splitThreshold3wayNs_ = std::atof(e) * 0.001;  // ps → ns
  }

  // Tier-aware cluster split
  if (const char* e = std::getenv("CTS_ENABLE_TIER_SPLIT"))
    enableTierSplit_ = (std::atoi(e) != 0);
  logger_->info(CTS, 376, "Tier split: {}", enableTierSplit_ ? "ON" : "OFF");

  // TARP clustering
  if (const char* e = std::getenv("CTS_ENABLE_TARP"))
    enableTarp_ = (std::atoi(e) != 0);
  if (const char* e = std::getenv("CTS_TARP_ALPHA"))
    tarpAlpha_ = std::stod(e);
  if (const char* e = std::getenv("CTS_ENABLE_PER_TIER_CTS"))
    enablePerTierCts_ = (std::atoi(e) != 0);
  logger_->info(CTS, 840, "TARP: {} (alpha={:.2f}), per-tier: {}",
                enableTarp_ ? "ON" : "OFF", tarpAlpha_,
                enablePerTierCts_ ? "ON" : "OFF");

  // Per-tier independent CTS.
  // MUST be placed BEFORE initSinkRegion() to avoid double preSinkClustering().
  // Env vars for grouped delay, relay, GH-tree etc. are read inside
  // readSkewDeliveryEnvVars() called from runSingleTierTree().
  if (enablePerTierCts_ && cts3dDb_ != nullptr) {
    // Read delivery env vars here (before runSingleTierTree) so they are set
    readSkewDeliveryEnvVars();

    int tier0Count = 0, tier1Count = 0;
    clock_.forEachSink([&](const ClockInst& sink) {
      odb::dbInst* dbInst = sink.getDbInst();
      if (dbInst == nullptr) return;
      int t = cts3dDb_->getInstTier(dbInst);
      if (t == 0) tier0Count++;
      else if (t == 1) tier1Count++;
    });
    logger_->info(CTS, 850, "Per-tier split: tier0={} FFs, tier1={} FFs",
                  tier0Count, tier1Count);

    // Create ONE shared root buffer for both tiers
    createSharedRoot();

    // Build per-tier sub-trees (L1+ branches attach to shared root)
    if (tier0Count > 0) runSingleTierTree(0);
    if (tier1Count > 0) runSingleTierTree(1);
    return;  // Skip unified tree path
  }

  auto prof_t0 = std::chrono::high_resolution_clock::now();
  initSinkRegion();
  auto prof_t1 = std::chrono::high_resolution_clock::now();
  double prof_clustering_s = std::chrono::duration<double>(prof_t1 - prof_t0).count();
  logger_->info(CTS, 911, "Profile: initSinkRegion (TARP+clustering) {:.1f}s", prof_clustering_s);

  // Fix G: Macro flat bypass (same logic as runSingleTierTree path)
  {
    int macroFlatMax = 50;
    if (const char* e = std::getenv("CTS_MACRO_FLAT_MAX_SINKS")) {
      macroFlatMax = std::atoi(e);
    }
    const int nSinks = clock_.getNumSinks();
    if (type_ == TreeType::MacroTree && macroFlatMax > 0 && nSinks <= macroFlatMax) {
      bool allBlock = true;
      clock_.forEachSink([&](const ClockInst& sink) {
        odb::dbInst* dbInst = sink.getDbInst();
        if (!isMacroBlockInst(dbInst)) {
          allBlock = false;
        }
      });
      if (allBlock) {
        logger_->info(CTS, 695,
            "Fix G: Macro flat bypass — {} sinks (all block), threshold={}, "
            "skipping H-tree → flat topology",
            nSinks, macroFlatMax);
        createSingleBufferClockNet();
        treeBufLevels_++;
        return;
      }
    }
  }

  // V58 Fix E: Macro tree depth clamp (fallback for non-flat macro trees)
  int effectiveMaxDepth = clockTreeMaxDepth_;
  if (type_ == TreeType::MacroTree) {
    int macroMaxDepth = 2;
    if (const char* e = std::getenv("CTS_MACRO_TREE_MAX_DEPTH")) {
      macroMaxDepth = std::atoi(e);
    }
    if (macroMaxDepth > 0 && macroMaxDepth < effectiveMaxDepth) {
      logger_->info(CTS, 692,
          "V58 Fix E: Macro tree depth clamped {} -> {} (sinks={})",
          effectiveMaxDepth, macroMaxDepth, clock_.getNumSinks());
      effectiveMaxDepth = macroMaxDepth;
    }
  }

  auto prof_htree_t0 = std::chrono::high_resolution_clock::now();
  for (int level = 1; level <= effectiveMaxDepth; ++level) {
    const unsigned numSinksPerSubRegion
        = computeNumberOfSinksPerSubRegion(level);
    double regionWidth, regionHeight;
    computeSubRegionSize(level, regionWidth, regionHeight);

    if (isSubRegionTooSmall(regionWidth, regionHeight)) {
      if (options_->isFakeLutEntriesEnabled()) {
        const unsigned minIndex = 1;
        techChar_->createFakeEntries(minLengthSinkRegion_, minIndex);
        minLengthSinkRegion_ = 1;
      } else {
        logger_->info(
            CTS,
            31,
            " Stop criterion found. Min length of sink region is ({}).",
            minLengthSinkRegion_);
        break;
      }
    }

    computeLevelTopology(level, regionWidth, regionHeight);

    if (isNumberOfSinksTooSmall(numSinksPerSubRegion)) {
      logger_->info(CTS,
                    32,
                    " Stop criterion found. Max number of sinks is {}.",
                    options_->getMaxFanout() ? options_->getMaxFanout()
                                             : numMaxLeafSinks_);
      break;
    }
  }

  auto prof_htree_t1 = std::chrono::high_resolution_clock::now();
  double prof_htree_s = std::chrono::duration<double>(prof_htree_t1 - prof_htree_t0).count();
  logger_->info(CTS, 912, "Profile: H-tree level loop ({} levels) {:.1f}s",
                static_cast<int>(topologyForEachLevel_.size()), prof_htree_s);

  if (topologyForEachLevel_.empty()) {
    createSingleBufferClockNet();
    treeBufLevels_++;
    return;
  }

  clock_.setMaxLevel(topologyForEachLevel_.size());

  if (options_->getPlotSolution()
      || logger_->debugCheck(utl::CTS, "HTree", 2)) {
    plotSolution();
  }

  if (CtsObserver* observer = options_->getObserver()) {
    observer->initializeWithClock(this, clock_);
  }

  if (options_->getObstructionAware()) {
    legalize();
  }

  // Extracted to readSkewDeliveryEnvVars().
  // Called from per-tier block (above) or here for unified mode.
  readSkewDeliveryEnvVars();

  auto prof_subnets_t0 = std::chrono::high_resolution_clock::now();
  createClockSubNets();
  auto prof_subnets_t1 = std::chrono::high_resolution_clock::now();
  double prof_subnets_s = std::chrono::duration<double>(prof_subnets_t1 - prof_subnets_t0).count();
  logger_->info(CTS, 913, "Profile: createClockSubNets (leaf+TAP+delay) {:.1f}s", prof_subnets_s);

  // Total run() time
  auto prof_total_t1 = std::chrono::high_resolution_clock::now();
  double prof_total_s = std::chrono::duration<double>(prof_total_t1 - prof_t0).count();
  logger_->info(CTS, 914,
      "Profile: run() total {:.1f}s (clustering={:.1f} htree={:.1f} subnets={:.1f})",
      prof_total_s, prof_clustering_s, prof_htree_s, prof_subnets_s);

  // clang-format off
  debugPrint(logger_, CTS, "legalizer", 3, "Htree file {} has been generated",
             plotHTree());
  debugPrint(logger_, CTS, "legalizer", 3, "Run 'obsAwareCts.py cts.clk.buffer'"
	     "to produce cts.clk.buffer.png");
  // clang-format on
}

// Extracted from run() — reads all skew delivery env vars.
// Safe to call multiple times (idempotent).
void HTreeBuilder::readSkewDeliveryEnvVars()
{
  if (cts3dDb_ == nullptr || !cts3dDb_->hasSkewTargets()) return;

  const char* thr_env = std::getenv("CTS_DELAY_TARGET_THRESHOLD");
  delayTargetThreshold_ = (thr_env != nullptr) ? std::atof(thr_env) : 0.010;

    // Single delay-buffer delay (ns) and max chain length per leaf
    if (const char* e = std::getenv("CTS_LEAF_BUF_DELAY_NS"))
      singleBufDelay_ = std::atof(e);
    if (const char* e = std::getenv("CTS_MAX_LEAF_DELAY_BUFS"))
      maxLeafDelayBufs_ = std::atoi(e);

    // Target-aware sub-clustering parameters
    // NOTE: env var reading moved to before initSinkRegion() (V42-fix).
    // Only logging remains here.
    if (enableTargetSplit_) {
      logger_->info(CTS, 410,
                    "Target-aware sub-clustering enabled: "
                    "2way_threshold={:.1f}ps, 3way_threshold={:.1f}ps",
                    splitThreshold2wayNs_ * 1000.0,
                    splitThreshold3wayNs_ * 1000.0);
    }

    // Per-FF relay buffer parameters
    enablePerFfRelay_ = false;
    if (const char* e = std::getenv("CTS_ENABLE_PER_FF_RELAY"))
      enablePerFfRelay_ = (std::atoi(e) != 0);
    if (enablePerFfRelay_ && cts3dDb_ && cts3dDb_->hasSkewTargets()) {
      if (const char* e = std::getenv("CTS_PER_FF_BUF_DELAY_NS"))
        perFfBufDelay_ = std::atof(e);
      if (const char* e = std::getenv("CTS_PER_FF_MAX_RELAY"))
        perFfMaxRelay_ = std::atoi(e);
      if (const char* e = std::getenv("CTS_PER_FF_HOLD_GUARD_NS"))
        perFfHoldGuard_ = std::atof(e);
      if (const char* e = std::getenv("CTS_PER_FF_MIN_TARGET_NS"))
        perFfMinTarget_ = std::atof(e);
      {
        auto hb_t0 = std::chrono::high_resolution_clock::now();
        if (const char* csv = std::getenv("CTS_TIMING_GRAPH_CSV"))
          loadPerFfHoldBudgets(csv);
        if (const char* ioCsv = std::getenv("CTS_IO_TIMING_CSV"))
          loadIoHoldBudgets(ioCsv);
        auto hb_t1 = std::chrono::high_resolution_clock::now();
        double hb_s = std::chrono::duration<double>(hb_t1 - hb_t0).count();
        logger_->info(CTS, 917, "Profile: loadHoldBudgets (CSV parse) {:.1f}s", hb_s);
      }
      // Elmore wire RC parameters for x_useful positioning
      relayRwKOhmPerUm_ = 0.0;
      relayCwFfPerUm_   = 0.0;
      relayDbuPerUm_    = 2000.0;
      if (const char* e = std::getenv("CTS_RELAY_RW_PER_UM"))
        relayRwKOhmPerUm_ = std::atof(e);
      if (const char* e = std::getenv("CTS_RELAY_CW_PER_UM"))
        relayCwFfPerUm_ = std::atof(e);
      if (const char* e = std::getenv("CTS_DBU_PER_UM"))
        relayDbuPerUm_ = std::atof(e);
      const bool useElmore = (relayRwKOhmPerUm_ > 1e-12 && relayCwFfPerUm_ > 1e-12);
      logger_->info(CTS, 420,
                    "Per-FF relay: bufDelay={:.3f}ns maxRelay={} "
                    "holdGuard={:.3f}ns minTarget={:.3f}ns holdBudgets={} "
                    "Elmore={}(rw={:.4f}kOhm/um cw={:.4f}fF/um)",
                    perFfBufDelay_, perFfMaxRelay_, perFfHoldGuard_,
                    perFfMinTarget_,
                    static_cast<int>(perFfHoldBudget_.size()),
                    useElmore ? "ON" : "OFF",
                    relayRwKOhmPerUm_, relayCwFfPerUm_);
    } else if (enablePerFfRelay_) {
      logger_->warn(CTS, 421,
                    "CTS_ENABLE_PER_FF_RELAY=1 but no skew targets; "
                    "falling back to per-branch mode");
      enablePerFfRelay_ = false;
    }

    // Grouped Delay Chain parameters.
    // Groups FFs within each leaf cluster by quantized LP target delta.
    // Shared delay chain per group → fewer nets than per-FF relay → avoids GRT-0183.
    enableGroupedDelay_ = false;
    if (const char* e = std::getenv("CTS_ENABLE_GROUPED_DELAY"))
      enableGroupedDelay_ = (std::atoi(e) != 0);
    if (enableGroupedDelay_ && cts3dDb_ && cts3dDb_->hasSkewTargets()) {
      if (const char* e = std::getenv("CTS_GROUPED_DELAY_MAX_DEPTH"))
        groupedDelayMaxDepth_ = std::atoi(e);
      if (const char* e = std::getenv("CTS_GROUPED_DELAY_MIN_DELTA_NS"))
        groupedDelayMinDelta_ = std::atof(e);
      // Reuse per-FF relay parameters: perFfBufDelay_, perFfHoldGuard_, perFfHoldBudget_
      // Load hold budgets if not already loaded by per-FF relay block
      if (perFfHoldBudget_.empty()) {
        if (const char* e = std::getenv("CTS_PER_FF_BUF_DELAY_NS"))
          perFfBufDelay_ = std::atof(e);
        if (const char* e = std::getenv("CTS_PER_FF_HOLD_GUARD_NS"))
          perFfHoldGuard_ = std::atof(e);
        if (const char* e = std::getenv("CTS_PER_FF_MIN_TARGET_NS"))
          perFfMinTarget_ = std::atof(e);
        {
          auto hb2_t0 = std::chrono::high_resolution_clock::now();
          if (const char* csv = std::getenv("CTS_TIMING_GRAPH_CSV"))
            loadPerFfHoldBudgets(csv);
          if (const char* ioCsv = std::getenv("CTS_IO_TIMING_CSV"))
            loadIoHoldBudgets(ioCsv);
          auto hb2_t1 = std::chrono::high_resolution_clock::now();
          double hb2_s = std::chrono::duration<double>(hb2_t1 - hb2_t0).count();
          logger_->info(CTS, 918, "Profile: loadHoldBudgets (per-tier path) {:.1f}s", hb2_s);
        }
      }
      // Cluster-Uniform Depth mode — one chain per cluster, depth = cluster median.
      enableClusterUniform_ = false;
      if (const char* e = std::getenv("CTS_GROUPED_DELAY_CLUSTER_UNIFORM"))
        enableClusterUniform_ = (std::atoi(e) != 0);

      // Liberty-based buffer delay override (from V51-BUF).
      // charBuf_ driving its own input cap = cascaded tap scenario.
      // Overrides env var CTS_PER_FF_BUF_DELAY_NS with Liberty-measured value.
      {
        const double libertyDelay = techChar_->getCharBufDelay();
        const double libInputCap  = techChar_->getCharBufInputCap();
        if (libertyDelay > 0) {
          perFfBufDelay_ = libertyDelay;
          logger_->info(CTS, 860,
              "Liberty buf delay: {:.4f}ns (charBuf input_cap={:.4f}pF)",
              perFfBufDelay_, libInputCap);
        } else {
          logger_->warn(CTS, 861,
              "Liberty delay not available; using default perFfBufDelay={:.3f}ns",
              perFfBufDelay_);
        }
      }

      // GH-Tree Cascaded Delay Tap.
      // When enabled, buildGroupChains uses cascaded tap structure instead of
      // parallel chains. Same depth delivery, ~50% fewer buffers and nets.
      enableGhTree_ = false;
      if (const char* e = std::getenv("CTS_ENABLE_GH_TREE"))
        enableGhTree_ = (std::atoi(e) != 0);
      // Adaptive depth cost ratio (default 1.0 = 1 FF per empty tap to justify)
      if (const char* e = std::getenv("CTS_ADAPTIVE_DEPTH_COST_RATIO"))
        adaptiveDepthCostRatio_ = std::atof(e);

      // Log after all env vars are read so values are accurate
      logger_->info(CTS, 490,
                    "Grouped delay: maxDepth={} minDelta={:.3f}ns "
                    "bufDelay={:.3f}ns holdBudgets={} "
                    "clusterUniform={} cascadedTap={} adaptiveCost={:.1f}",
                    groupedDelayMaxDepth_, groupedDelayMinDelta_,
                    perFfBufDelay_,
                    static_cast<int>(perFfHoldBudget_.size()),
                    enableClusterUniform_ ? 1 : 0,
                    enableGhTree_ ? 1 : 0,
                    adaptiveDepthCostRatio_);
      if (enableGhTree_) {
        logger_->info(CTS, 570,
            "Cascaded tap + adaptive depth: costRatio={:.1f}, bufDelay={:.3f}ns",
            adaptiveDepthCostRatio_, perFfBufDelay_);
      }
    } else if (enableGroupedDelay_) {
      logger_->warn(CTS, 491,
                    "CTS_ENABLE_GROUPED_DELAY=1 but no skew targets; "
                    "falling back to balanced CTS");
      enableGroupedDelay_ = false;
    }

    // Intermediate delay buffer gate (V37 default on, V39 default off)
    enableMidDelayBufs_ = false;
    if (const char* e = std::getenv("CTS_ENABLE_MID_DELAY_BUFS"))
      enableMidDelayBufs_ = (std::atoi(e) != 0);

    // CKMeans target-aware branching weight.
    // When > 0, H-tree 2-way split (CKMeans) groups FFs with similar targets.
    ckmeansTargetBeta_ = 0.0;
    if (const char* e = std::getenv("CTS_CKMEANS_TARGET_BETA"))
      ckmeansTargetBeta_ = std::atof(e);
    if (ckmeansTargetBeta_ > 0.0) {
      logger_->info(CTS, 428,
          "CKMeans target-aware branching: beta={:.1f}",
          ckmeansTargetBeta_);
    }

    computeBranchDelayTargets();
}

bool HTreeBuilder::isNumberOfSinksTooSmall(unsigned numSinksPerSubRegion) const
{
  if (options_->getMaxFanout()) {
    return numSinksPerSubRegion < options_->getMaxFanout();
  }
  return numSinksPerSubRegion < numMaxLeafSinks_;
}

std::string HTreeBuilder::plotHTree()
{
  auto name = std::string("cts.") + clock_.getName() + ".buffer";
  std::ofstream file(name);

  plotBlockage(file, db_, wireSegmentUnit_);

  Point<double> topLevelBufferLoc = sinkRegion_.getCenter();

  for (int levelIdx = 0; levelIdx < topologyForEachLevel_.size(); ++levelIdx) {
    LevelTopology& topology = topologyForEachLevel_[levelIdx];

    // clang-format off
    topology.forEachBranchingPoint(
        [&](unsigned idx, Point<double> branchPoint) {
          unsigned parentIdx = topology.getBranchingPointParentIdx(idx);

          Point<double> parentPoint
              = (levelIdx == 0)
                    ? topLevelBufferLoc
                    : topologyForEachLevel_[levelIdx - 1].getBranchingPoint(
                          parentIdx);

          const std::vector<Point<double>>& sinks
              = topology.getBranchSinksLocations(idx);

          plotSinks(file, sinks);

          double x1 = parentPoint.getX();
          double y1 = parentPoint.getY();
          double x2 = branchPoint.getX();
          double y2 = branchPoint.getY();
          std::string name = "buffer";
          file << levelIdx << " " << x1 << " " << y1 << " " << x2 << " " << y2;
          file << " " << name << '\n';
        });
    // clang-format on
  }

  LevelTopology& leafTopology = topologyForEachLevel_.back();
  unsigned numSinks = 0;
  leafTopology.forEachBranchingPoint(
      [&](unsigned idx, Point<double> branchPoint) {
        double px = branchPoint.getX();
        double py = branchPoint.getY();

        const std::vector<Point<double>>& sinkLocs
            = leafTopology.getBranchSinksLocations(idx);

        for (const Point<double>& loc : sinkLocs) {
          auto name2 = mapLocationToSink_[loc]->getName();

          file << numSinks << " " << loc.getX() << " " << loc.getY();
          file << " " << px << " " << py << " leafbuffer " << name2;
          file << " z=" << wireSegmentUnit_ << '\n';
          ++numSinks;
        }
      });
  file.close();
  return name;
}

unsigned HTreeBuilder::computeNumberOfSinksPerSubRegion(
    const unsigned level) const
{
  unsigned min_clustering_sinks = (type_ == TreeType::MacroTree)
                                      ? min_clustering_macro_sinks_
                                      : min_clustering_sinks_;
  unsigned totalNumSinks = 0;
  if (clock_.getNumSinks() > min_clustering_sinks
      && options_->getSinkClustering()) {
    totalNumSinks = topLevelSinksClustered_.size();
  } else {
    totalNumSinks = clock_.getNumSinks();
  }
  const unsigned numRoots = std::pow(2, level);
  const double numSinksPerRoot = (double) totalNumSinks / numRoots;
  return std::ceil(numSinksPerRoot);
}

void HTreeBuilder::computeSubRegionSize(const unsigned level,
                                        double& width,
                                        double& height) const
{
  unsigned gridSizeX = 0;
  unsigned gridSizeY = 0;
  if (isVertical(1)) {
    gridSizeY = computeGridSizeX(level);
    gridSizeX = computeGridSizeY(level);
  } else {
    gridSizeX = computeGridSizeX(level);
    gridSizeY = computeGridSizeY(level);
  }
  width = sinkRegion_.getWidth() / gridSizeX;
  height = sinkRegion_.getHeight() / gridSizeY;
}

void HTreeBuilder::computeLevelTopology(const unsigned level,
                                        const double width,
                                        const double height)
{
  const unsigned numSinksPerSubRegion = computeNumberOfSinksPerSubRegion(level);
  logger_->report(" Level {}", level);
  logger_->report("    Direction: {}",
                  (isVertical(level)) ? ("Vertical") : ("Horizontal"));
  logger_->report("    Sinks per sub-region: {}", numSinksPerSubRegion);
  logger_->report("    Sub-region size: {:.4f} X {:.4f}", width, height);

  const unsigned minLength = minLengthSinkRegion_;
  const unsigned clampedMinLength = std::max(minLength, 1u);

  unsigned segmentLength
      = std::round(width / (double) clampedMinLength) * minLength / 2;

  if (isVertical(level)) {
    segmentLength
        = std::round(height / (double) clampedMinLength) * minLength / 2;
  }
  segmentLength = std::max<unsigned>(segmentLength, 1);

  LevelTopology topology(segmentLength);

  logger_->info(CTS, 34, "    Segment length (rounded): {}.", segmentLength);

  const int vertexBufferLength
      = options_->getVertexBufferDistance() / (techChar_->getLengthUnit() * 2);
  int remainingLength
      = options_->getBufferDistance() / (techChar_->getLengthUnit());
  int currWl = 0;
  unsigned inputCap = minInputCap_;
  unsigned inputSlew = 1;
  if (level > 1) {
    const LevelTopology& previousLevel = topologyForEachLevel_[level - 2];
    inputCap = previousLevel.getOutputCap();
    inputSlew = previousLevel.getOutputSlew();
    remainingLength = previousLevel.getRemainingLength();
    currWl = previousLevel.getCurrWl();
  }

  const unsigned kSlewThreshold = options_->getMaxSlew();
  const unsigned kInitTolerance = 1;

  int wirelengthThreshold;
  // If max wirelength is 0, set it as slew threshold  * maximum topology
  // wirelength. This will behave as if there was no max wirelength threshold.
  if (!options_->getMaxWl()) {
    wirelengthThreshold = kSlewThreshold * techChar_->getMaxSegmentLength();
  } else {
    wirelengthThreshold = options_->getMaxWl() / options_->getWireSegmentUnit();
  }

  debugPrint(
      logger_, CTS, "tech char", 1, "slew threshold = {}", kSlewThreshold);
  debugPrint(logger_,
             CTS,
             "tech char",
             1,
             "wirelength threshold = {}",
             wirelengthThreshold);
  unsigned length = 0;
  for (int charSegLength = techChar_->getMaxSegmentLength(); charSegLength >= 1;
       --charSegLength) {
    const unsigned numWires = (segmentLength - length) / charSegLength;

    if (numWires >= 1) {
      for (int wireCount = 0; wireCount < numWires; ++wireCount) {
        debugPrint(logger_, CTS, "tech char", 1, "curr wl = {}", currWl);
        unsigned outCap = 0, outSlew = 0;
        unsigned key = 0;
        if (options_->isSimpleSegmentEnabled()) {
          remainingLength -= charSegLength;

          if (segmentLength >= vertexBufferLength && (wireCount + 1 >= numWires)
              && options_->isVertexBuffersEnabled()) {
            remainingLength = 0;
            key = computeMinDelaySegment(charSegLength,
                                         inputSlew,
                                         inputCap,
                                         kSlewThreshold,
                                         kInitTolerance,
                                         outSlew,
                                         outCap,
                                         true,
                                         remainingLength);
            remainingLength
                += options_->getBufferDistance() / (techChar_->getLengthUnit());
          } else {
            if (remainingLength <= 0) {
              key = computeMinDelaySegment(charSegLength,
                                           inputSlew,
                                           inputCap,
                                           kSlewThreshold,
                                           kInitTolerance,
                                           outSlew,
                                           outCap,
                                           true,
                                           remainingLength);
              remainingLength += options_->getBufferDistance()
                                 / (techChar_->getLengthUnit());
            } else {
              key = computeMinDelaySegment(charSegLength,
                                           inputSlew,
                                           inputCap,
                                           kSlewThreshold,
                                           kInitTolerance,
                                           outSlew,
                                           outCap,
                                           false,
                                           remainingLength);
            }
          }
        } else {
          key = computeMinDelaySegment(charSegLength,
                                       inputSlew,
                                       inputCap,
                                       kSlewThreshold,
                                       wirelengthThreshold,
                                       kInitTolerance,
                                       outSlew,
                                       outCap,
                                       currWl);
        }

        if (key == std::numeric_limits<unsigned>::max()) {
          // No tech char entry found.
          continue;
        }

        length += charSegLength;
        techChar_->reportSegment(key);

        inputCap = std::max(outCap, minInputCap_);
        inputSlew = outSlew;
        topology.addWireSegment(key);
        topology.setRemainingLength(remainingLength);
      }

      if (length == segmentLength) {
        break;
      }
    }
  }

  topology.setOutputSlew(inputSlew);
  topology.setOutputCap(inputCap);
  topology.setCurrWl(currWl);

  computeBranchingPoints(level, topology);
  topologyForEachLevel_.push_back(topology);
}

unsigned HTreeBuilder::computeMinDelaySegment(const unsigned length) const
{
  unsigned minKey = std::numeric_limits<unsigned>::max();
  unsigned minDelay = std::numeric_limits<unsigned>::max();

  techChar_->forEachWireSegment(
      length, 1, 1, [&](unsigned key, const WireSegment& seg) {
        if (!seg.isBuffered()) {
          return;
        }
        if (seg.getDelay() < minDelay) {
          minKey = key;
          minDelay = seg.getDelay();
        }
      });

  return minKey;
}

unsigned HTreeBuilder::computeMinDelaySegment(const unsigned length,
                                              const unsigned inputSlew,
                                              const unsigned inputCap,
                                              const unsigned slewThreshold,
                                              const int wirelengthThreshold,
                                              const unsigned tolerance,
                                              unsigned& outputSlew,
                                              unsigned& outputCap,
                                              int& currWl) const
{
  unsigned minKey = std::numeric_limits<unsigned>::max();
  unsigned minDelay = std::numeric_limits<unsigned>::max();
  unsigned minBufKey = std::numeric_limits<unsigned>::max();
  unsigned minBufDelay = std::numeric_limits<unsigned>::max();

  for (int load = 1; load <= techChar_->getMaxCapacitance(); ++load) {
    for (int outSlew = 1; outSlew <= techChar_->getMaxSlew(); ++outSlew) {
      techChar_->forEachWireSegment(
          length, load, outSlew, [&](unsigned key, const WireSegment& seg) {
            if (std::abs((int) seg.getInputCap() - (int) inputCap) > tolerance
                || std::abs((int) seg.getInputSlew() - (int) inputSlew)
                       > tolerance) {
              return;
            }

            if (seg.isBuffered()
                && (currWl + seg.getWl2FirstBuffer() > wirelengthThreshold)) {
              return;
            }

            if (!seg.isBuffered()
                && (currWl + seg.getWl2FirstBuffer()
                    > wirelengthThreshold - techChar_->getMinSegmentLength())) {
              return;
            }

            if (seg.getDelay() < minDelay) {
              minDelay = seg.getDelay();
              minKey = key;
            }

            if (seg.isBuffered() && seg.getDelay() < minBufDelay) {
              minBufDelay = seg.getDelay();
              minBufKey = key;
            }
          });
    }
  }

  const unsigned MAX_TOLERANCE = 10;
  if (inputSlew >= slewThreshold) {
    if (minBufKey < std::numeric_limits<unsigned>::max()) {
      const WireSegment& bestBufSegment = techChar_->getWireSegment(minBufKey);
      outputSlew = bestBufSegment.getOutputSlew();
      outputCap = bestBufSegment.getLoad();
      currWl = bestBufSegment.getLastWl();
      return minBufKey;
    }
    if (tolerance < MAX_TOLERANCE) {
      // Increasing tolerance
      return computeMinDelaySegment(length,
                                    inputSlew,
                                    inputCap,
                                    slewThreshold,
                                    wirelengthThreshold,
                                    tolerance + 1,
                                    outputSlew,
                                    outputCap,
                                    currWl);
    }
  }

  if (minKey == std::numeric_limits<unsigned>::max()) {
    if (tolerance >= MAX_TOLERANCE) {
      return minKey;
    }
    // Increasing tolerance
    return computeMinDelaySegment(length,
                                  inputSlew,
                                  inputCap,
                                  slewThreshold,
                                  wirelengthThreshold,
                                  tolerance + 1,
                                  outputSlew,
                                  outputCap,
                                  currWl);
  }

  const WireSegment& bestSegment = techChar_->getWireSegment(minKey);
  if (bestSegment.isBuffered()) {
    outputSlew = bestSegment.getOutputSlew();
    currWl = bestSegment.getLastWl();
  } else {
    outputSlew
        = std::max((unsigned) bestSegment.getOutputSlew(), inputSlew + 1);
    currWl += bestSegment.getLastWl();
  }
  outputCap = bestSegment.getLoad();

  return minKey;
}

unsigned HTreeBuilder::computeMinDelaySegment(const unsigned length,
                                              const unsigned inputSlew,
                                              const unsigned inputCap,
                                              const unsigned slewThreshold,
                                              const unsigned tolerance,
                                              unsigned& outputSlew,
                                              unsigned& outputCap,
                                              const bool forceBuffer,
                                              const int expectedLength) const
{
  unsigned minKey = std::numeric_limits<unsigned>::max();
  unsigned minDelay = std::numeric_limits<unsigned>::max();
  unsigned minBufKey = std::numeric_limits<unsigned>::max();
  unsigned minBufDelay = std::numeric_limits<unsigned>::max();
  unsigned minBufKeyFallback = std::numeric_limits<unsigned>::max();
  unsigned minDelayFallback = std::numeric_limits<unsigned>::max();

  for (int load = 1; load <= techChar_->getMaxCapacitance(); ++load) {
    for (int outSlew = 1; outSlew <= techChar_->getMaxSlew(); ++outSlew) {
      techChar_->forEachWireSegment(
          length, load, outSlew, [&](unsigned key, const WireSegment& seg) {
            // Same as the other functions, however, forces a segment
            // to have a buffer in a specific location.
            const unsigned normalLength = length;
            if (!seg.isBuffered() && seg.getDelay() < minDelay) {
              minDelay = seg.getDelay();
              minKey = key;
            }
            if (seg.isBuffered() && seg.getDelay() < minBufDelay
                && seg.getNumBuffers() == 1) {
              // If buffer is in the range of 10% of the expected location, save
              // its key.
              if (seg.getBufferLocation(0)
                      > ((((double) normalLength + (double) expectedLength)
                          / (double) normalLength)
                         * 0.9)
                  && seg.getBufferLocation(0)
                         < ((((double) normalLength + (double) expectedLength)
                             / (double) normalLength)
                            * 1.1)) {
                minBufDelay = seg.getDelay();
                minBufKey = key;
              }
              if (seg.getDelay() < minDelayFallback) {
                minDelayFallback = seg.getDelay();
                minBufKeyFallback = key;
              }
            }
          });
    }
  }

  if (forceBuffer) {
    if (minBufKey != std::numeric_limits<unsigned>::max()) {
      return minBufKey;
    }

    if (minBufKeyFallback != std::numeric_limits<unsigned>::max()) {
      return minBufKeyFallback;
    }
  }

  return minKey;
}

void HTreeBuilder::computeBranchingPoints(const unsigned level,
                                          LevelTopology& topology)
{
  if (level == 1) {
    const Point<double> clockRoot(sinkRegion_.getCenter());
    Point<double> low(clockRoot);
    Point<double> high(clockRoot);
    if (isHorizontal(level)) {
      low.setX(low.getX() - topology.getLength());
      high.setX(high.getX() + topology.getLength());
    } else {
      low.setY(low.getY() - topology.getLength());
      high.setY(high.getY() + topology.getLength());
    }
    const unsigned branchPtIdx1
        = topology.addBranchingPoint(low, LevelTopology::NO_PARENT);
    const unsigned branchPtIdx2
        = topology.addBranchingPoint(high, LevelTopology::NO_PARENT);

    refineBranchingPointsWithClustering(topology,
                                        level,
                                        branchPtIdx1,
                                        branchPtIdx2,
                                        clockRoot,
                                        topLevelSinksClustered_);
    return;
  }

  LevelTopology& parentTopology = topologyForEachLevel_[level - 2];
  parentTopology.forEachBranchingPoint(
      [&](unsigned idx, Point<double> clockRoot) {
        Point<double> low(clockRoot);
        Point<double> high(clockRoot);
        if (isHorizontal(level)) {
          low.setX(low.getX() - topology.getLength());
          high.setX(high.getX() + topology.getLength());
        } else {
          low.setY(low.getY() - topology.getLength());
          high.setY(high.getY() + topology.getLength());
        }
        const unsigned branchPtIdx1 = topology.addBranchingPoint(low, idx);
        const unsigned branchPtIdx2 = topology.addBranchingPoint(high, idx);

        std::vector<std::pair<float, float>> sinks;
        computeBranchSinks(parentTopology, idx, sinks);
        refineBranchingPointsWithClustering(
            topology, level, branchPtIdx1, branchPtIdx2, clockRoot, sinks);
      });
}

// Per-tier independent H-tree.
// Create a single shared root buffer at the center of ALL sinks (both tiers).
// Both per-tier sub-trees attach their L1 branches to this shared root subnet.
// This ensures tier-to-tier clock skew is controlled at the root level.
void HTreeBuilder::createSharedRoot()
{
  wireSegmentUnit_ = techChar_->getLengthUnit();

  // Compute center from ALL sinks (not filtered by tier)
  std::vector<std::pair<float, float>> allSinks;
  std::vector<const ClockInst*> allInsts;
  initTopLevelSinks(allSinks, allInsts);

  float minX = 1e30f, maxX = -1e30f, minY = 1e30f, maxY = -1e30f;
  for (auto& [x, y] : allSinks) {
    minX = std::min(minX, x); maxX = std::max(maxX, x);
    minY = std::min(minY, y); maxY = std::max(maxY, y);
  }
  Point<double> center((minX + maxX) / 2.0, (minY + maxY) / 2.0);

  // Root tier: overall dominant tier from all sinks.
  // V58_CL2: CTS_FORCE_ROOT_TIER overrides dominant-tier selection.
  // Pin3D uses bottom root buffer (same tier as majority of macros).
  // Our getDominantTierFromClock() picks upper (98% upper FFs in swerv),
  // creating a fundamentally different tree from Pin3D even with per-tier OFF.
  int root_tier = cts3dDb_->getDominantTierFromClock(clock_);
  if (const char* e = std::getenv("CTS_FORCE_ROOT_TIER")) {
    int forced = std::atoi(e);
    if (forced >= 0 && forced != root_tier) {
      logger_->info(CTS, 700,
          "V58_CL2: Root tier forced {} -> {} (CTS_FORCE_ROOT_TIER)",
          root_tier, forced);
      root_tier = forced;
    }
  }
  sharedRootTier_ = root_tier;
  const std::string root_buffer
      = cts3dDb_->getBufferForTier(options_->getRootBuffer(), root_tier);

  Point<double> legalCenter = legalizeOneBuffer(center, root_buffer);
  const int centerX = legalCenter.getX() * wireSegmentUnit_;
  const int centerY = legalCenter.getY() * wireSegmentUnit_;

  // Create shared root buffer and subnet (no tier prefix — shared by both tiers)
  sharedRootBuffer_
      = &clock_.addClockBuffer("clkbuf_0", root_buffer, centerX, centerY);
  // Only set topBufferName_ if not already set by forkRegisterClockNetwork().
  // RegisterTree builders have topBufferName_ = "clkbuf_regs_0_..." (separator buffer
  // that connects macro clock net to register clock net). Overwriting it with the
  // shared root name breaks clock connectivity for register trees in multi-clock designs.
  if (topBufferName_.empty()) {
    topBufferName_ = sharedRootBuffer_->getName();
  }
  addTreeLevelBuffer(sharedRootBuffer_);

  sharedRootSubNet_ = &clock_.addSubNet("clknet_0");
  sharedRootSubNet_->setTreeLevel(0);
  sharedRootSubNet_->addInst(*sharedRootBuffer_);
  treeBufLevels_++;

  sharedRootLocation_ = legalCenter;

  logger_->info(CTS, 856,
      "Shared root: tier={}, buffer={}, sinks={}, center=({:.2f}, {:.2f}), "
      "rootBufName='{}', topBufferName='{}', clockName='{}'",
      root_tier, root_buffer, allSinks.size(),
      legalCenter.getX(), legalCenter.getY(),
      sharedRootBuffer_->getName(), topBufferName_, clock_.getName());
}

// Runs the complete H-tree pipeline for a single tier:
//   1. Filter sinks to this tier
//   2. Clear per-tree state (mapLocationToSink_, topLevelSinksClustered_, topology)
//   3. Run clustering → CKMeans → wire segments → createClockSubNets
// When sharedRootBuffer_ is set, L1 branches attach to the shared root
// instead of creating a tier-specific root buffer.
void HTreeBuilder::runSingleTierTree(int tier)
{
  logger_->info(CTS, 851, "Building independent H-tree for tier {}", tier);

  // --- Step 0: Set tier naming prefix to avoid buffer name collisions ---
  // Tier 0 uses no prefix (backward compatible names: clkbuf_0, clkbuf_leaf_0, etc.)
  // Tier 1 uses "t1_" prefix: t1_clkbuf_0, t1_clkbuf_leaf_0, etc.
  tierPrefix_ = (tier == 0) ? "" : "t" + std::to_string(tier) + "_";

  // --- Step 1: Clear per-tree state ---
  mapLocationToSink_.clear();
  topLevelSinksClustered_.clear();
  topologyForEachLevel_.clear();
  clusterBufTarget_.clear();
  branchDelayTargets_.clear();
  perLevelGlobalMean_.clear();

  // wireSegmentUnit_ comes from TechChar (charBuf height × 10 / 2).
  // charBuf is selected from buf_list (highest max-cap buffer).
  // Do NOT override per-tier: in heterogeneous 3D with upper-only buf_list,
  // using bottom tier buffer height (5x larger for NG45) destroys H-tree grid.
  // Removed per-tier wireSegmentUnit recalculation (V57 hotfix, 2026-04-03).

  // --- Step 2: Filter sinks for this tier and cluster ---
  double clusterDiameter = (type_ == TreeType::MacroTree)
                               ? options_->getMacroMaxDiameter()
                               : options_->getMaxDiameter();
  unsigned clusterSize = (type_ == TreeType::MacroTree)
                             ? options_->getMacroSinkClusteringSize()
                             : options_->getSinkClusteringSize();
  const int dbUnits = options_->getDbUnits();
  const float maxDiameter = (clusterDiameter * dbUnits) / wireSegmentUnit_;

  std::vector<std::pair<float, float>> topLevelSinks;
  std::vector<const ClockInst*> sinkInsts;
  initTopLevelSinksForTier(tier, topLevelSinks, sinkInsts);
  logger_->info(CTS, 852, "Tier {}: {} sinks after filtering", tier,
                topLevelSinks.size());

  if (topLevelSinks.empty()) {
    logger_->info(CTS, 853, "Tier {}: no sinks, skipping tree", tier);
    return;
  }

  // --- Step 4: Pre-sink clustering (Hilbert+Greedy or TARP) ---
  preSinkClustering(topLevelSinks, sinkInsts, maxDiameter, clusterSize);

  // Compute sink region from clustered sinks
  unsigned min_clustering_sinks = (type_ == TreeType::MacroTree)
                                      ? min_clustering_macro_sinks_
                                      : min_clustering_sinks_;
  if (topLevelSinks.size() <= min_clustering_sinks
      || !(options_->getSinkClustering())) {
    Box<int> sinkRegionDbu = clock_.computeSinkRegion();
    sinkRegion_ = sinkRegionDbu.normalize(1.0 / wireSegmentUnit_);
  } else {
    if (topLevelSinksClustered_.size() > 400
        && options_->getSinkClusteringLevels() > 0) {
      std::vector<std::pair<float, float>> secondLevelLocs;
      std::vector<const ClockInst*> secondLevelInsts;
      initSecondLevelSinks(secondLevelLocs, secondLevelInsts);
      preSinkClustering(secondLevelLocs, secondLevelInsts,
                        maxDiameter * 4,
                        std::ceil(std::sqrt(clusterSize)),
                        true);
    }
    sinkRegion_ = clock_.computeSinkRegionClustered(topLevelSinksClustered_);
  }
  logger_->info(CTS, 854, "Tier {} sink region: {}", tier, sinkRegion_);

  // --- Step 4b: Macro flat bypass (Fix G) ---
  // Pin3D: bottom-only CTS → macro net has ~7 sinks (<15) → single-buffer flat
  //   topology → macro latency ~563ps, balanced with registers (~442ps).
  // V58 per-tier CTS: all 21 SRAMs on clk net → H-tree builds → 1,095ps latency
  //   → 633ps gap vs registers → catastrophic hold violations.
  // Fix G: when MacroTree AND all sinks are block (SRAM) AND sink_count ≤ threshold,
  //   skip H-tree entirely → createSingleBufferClockNet (flat fan-out from root).
  //   This replicates Pin3D's behavior: root → direct wire → each macro.
  // CTS_MACRO_FLAT_MAX_SINKS env override (0=disable, default 50).
  {
    int macroFlatMax = 50;
    if (const char* e = std::getenv("CTS_MACRO_FLAT_MAX_SINKS")) {
      macroFlatMax = std::atoi(e);
    }
    const int nSinks = static_cast<int>(topLevelSinks.size());
    if (type_ == TreeType::MacroTree && macroFlatMax > 0 && nSinks <= macroFlatMax) {
      // Verify all sinks are block (SRAM/macro) instances
      bool allBlock = true;
      for (const auto* sinkInst : sinkInsts) {
        if (!sinkInst || !isMacroBlockInst(sinkInst->getDbInst())) {
          allBlock = false;
          break;
        }
      }
      if (allBlock) {
        logger_->info(CTS, 693,
            "Fix G: Macro flat bypass — {} sinks (all block), threshold={}, "
            "skipping H-tree → flat topology (Pin3D-equivalent)",
            nSinks, macroFlatMax);

        if (sharedRootSubNet_ != nullptr) {
          // Per-tier mode: shared root already exists (clkbuf_0_clk).
          // V58_CL Fix: Cross-tier macro sub-root buffer.
          // When macro sinks are on a different tier from the shared root,
          // connecting them directly creates long cross-tier wires with huge
          // latency spread (swerv DCCM: 475-1163ps from unbalanced wire RC).
          // Fix: create a per-tier sub-root buffer on the macro's tier,
          // connect it to the shared root, then connect macros to the sub-root.
          // Same-tier macros still connect directly (no change).
          const bool crossTier = (tier != sharedRootTier_
                                  && sharedRootTier_ >= 0);
          ClockSubNet* macroSubNet = sharedRootSubNet_;
          if (crossTier) {
            // Compute centroid of macro sinks on this tier
            float cx = 0, cy = 0;
            int cnt = 0;
            for (const ClockInst* inst : sinkInsts) {
              if (inst != nullptr) {
                cx += (double) inst->getX() / wireSegmentUnit_;
                cy += (double) inst->getY() / wireSegmentUnit_;
                cnt++;
              }
            }
            if (cnt > 0) { cx /= cnt; cy /= cnt; }
            // Create sub-root buffer on the macro's tier
            const std::string subRootMaster = cts3dDb_->getBufferForTier(
                options_->getRootBuffer(), tier);
            Point<double> subRootLoc(cx, cy);
            Point<double> legalSubRoot = legalizeOneBuffer(
                subRootLoc, subRootMaster);
            const int sx = legalSubRoot.getX() * wireSegmentUnit_;
            const int sy = legalSubRoot.getY() * wireSegmentUnit_;
            const std::string subRootName
                = tierPrefix_ + "clkbuf_macro_subroot_"
                  + clock_.getName();
            ClockInst& subRootBuf = clock_.addClockBuffer(
                subRootName, subRootMaster, sx, sy);
            if (cts3dDb_ != nullptr) {
              cts3dDb_->setClockInstTier(subRootBuf, tier);
            }
            addTreeLevelBuffer(&subRootBuf);
            // Connect sub-root to shared root subnet
            sharedRootSubNet_->addInst(subRootBuf);
            // Create sub-net for macro sinks
            const std::string subNetName
                = tierPrefix_ + "clknet_macro_subroot_"
                  + clock_.getName();
            ClockSubNet& subNet = clock_.addSubNet(subNetName);
            subNet.addInst(subRootBuf);
            subNet.setLeafLevel(true);
            macroSubNet = &subNet;
            treeBufLevels_++;
            logger_->info(CTS, 696,
                "Fix G: cross-tier macro sub-root '{}' on tier {} "
                "at ({:.2f}, {:.2f}) for {} macro sinks "
                "(shared root on tier {})",
                subRootName, tier, legalSubRoot.getX(),
                legalSubRoot.getY(), nSinks, sharedRootTier_);
          } else {
            // Same tier: connect directly to shared root (original behavior)
            if (!sharedRootSubNet_->isLeafLevel()) {
              sharedRootSubNet_->setLeafLevel(true);
              treeBufLevels_++;
            }
          }
          for (const ClockInst* inst : sinkInsts) {
            if (inst != nullptr) {
              macroSubNet->addInst(*const_cast<ClockInst*>(inst));
            }
          }
          logger_->info(CTS, 694,
              "Fix G: tier {} added {} macro sinks {} '{}'",
              tier, nSinks,
              crossTier ? "via sub-root to" : "directly to shared root",
              crossTier ? (tierPrefix_ + "clkbuf_macro_subroot_"
                           + clock_.getName())
                        : sharedRootBuffer_->getName());
        } else {
          // Non-per-tier mode: use standard single-buffer flat topology.
          createSingleBufferClockNet();
          treeBufLevels_++;
        }
        return;
      }
    }
  }

  // --- Step 5: H-tree level topology (CKMeans bisection) ---
  // V58 Fix E: Macro tree depth clamp (fallback for non-flat macro trees).
  int effectiveMaxDepth = clockTreeMaxDepth_;
  if (type_ == TreeType::MacroTree) {
    int macroMaxDepth = 2;
    if (const char* e = std::getenv("CTS_MACRO_TREE_MAX_DEPTH")) {
      macroMaxDepth = std::atoi(e);
    }
    if (macroMaxDepth > 0 && macroMaxDepth < effectiveMaxDepth) {
      logger_->info(CTS, 691,
          "V58 Fix E: Macro tree depth clamped {} -> {} (sinks={})",
          effectiveMaxDepth, macroMaxDepth,
          static_cast<int>(topLevelSinksClustered_.size()));
      effectiveMaxDepth = macroMaxDepth;
    }
  }
  // V58_CL: Per-tier register tree depth clamp.
  // Per-tier CTS creates deeper trees than single-tree Pin3D (path 14 vs 8)
  // because each tier builds an independent H-tree from the shared root.
  // Higher base latency → hold violations even for FFs with lp_target=0.
  // CTS_PER_TIER_REG_MAX_DEPTH (default 0=disabled) caps H-tree levels
  // for register trees in per-tier mode. Pin3D typically builds 7-8 levels.
  if (enablePerTierCts_ && type_ != TreeType::MacroTree) {
    int perTierRegMaxDepth = 0;
    if (const char* e = std::getenv("CTS_PER_TIER_REG_MAX_DEPTH")) {
      perTierRegMaxDepth = std::atoi(e);
    }
    if (perTierRegMaxDepth > 0
        && perTierRegMaxDepth < effectiveMaxDepth) {
      logger_->info(CTS, 697,
          "V58_CL: Per-tier register tree depth clamped {} -> {} "
          "(tier={}, sinks={})",
          effectiveMaxDepth, perTierRegMaxDepth, tier,
          static_cast<int>(topLevelSinksClustered_.size()));
      effectiveMaxDepth = perTierRegMaxDepth;
    }
  }
  for (int level = 1; level <= effectiveMaxDepth; ++level) {
    const unsigned numSinksPerSubRegion
        = computeNumberOfSinksPerSubRegion(level);
    double regionWidth, regionHeight;
    computeSubRegionSize(level, regionWidth, regionHeight);

    if (isSubRegionTooSmall(regionWidth, regionHeight)) {
      if (options_->isFakeLutEntriesEnabled()) {
        const unsigned minIndex = 1;
        techChar_->createFakeEntries(minLengthSinkRegion_, minIndex);
        minLengthSinkRegion_ = 1;
      } else {
        break;
      }
    }
    computeLevelTopology(level, regionWidth, regionHeight);
    if (isNumberOfSinksTooSmall(numSinksPerSubRegion)) {
      break;
    }
  }

  if (topologyForEachLevel_.empty()) {
    // Too few sinks for H-tree — create single buffer for this tier
    logger_->info(CTS, 855, "Tier {}: too few clusters for H-tree, "
                  "creating single buffer", tier);
    createSingleBufferClockNet();
    treeBufLevels_++;
    return;
  }

  clock_.setMaxLevel(topologyForEachLevel_.size());

  if (options_->getObstructionAware()) {
    legalize();
  }

  // --- Step 6: Compute branch delay targets (for grouped delay/TAP) ---
  if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()) {
    computeBranchDelayTargets();
  }

  // --- Step 7: Create clock sub-nets (buffers, wires, grouped delay) ---
  createClockSubNets();
}

void HTreeBuilder::initTopLevelSinks(
    std::vector<std::pair<float, float>>& sinkLocations,
    std::vector<const ClockInst*>& sinkInsts)
{
  sinkLocations.clear();
  clock_.forEachSink([&](const ClockInst& sink) {
    // FloatFix: (float) → (double) to prevent map key mismatch → CTS-0080
    sinkLocations.emplace_back((double) sink.getX() / wireSegmentUnit_,
                               (double) sink.getY() / wireSegmentUnit_);
    sinkInsts.emplace_back(&sink);
  });
}

// Filter sinks to a specific tier for per-tier CTS.
// Only includes FFs whose ODB instance is on the requested tier
// (determined by Cts3DDatabase::getInstTier).
void HTreeBuilder::initTopLevelSinksForTier(
    int tier,
    std::vector<std::pair<float, float>>& sinkLocations,
    std::vector<const ClockInst*>& sinkInsts)
{
  sinkLocations.clear();
  sinkInsts.clear();
  clock_.forEachSink([&](const ClockInst& sink) {
    if (cts3dDb_ == nullptr) return;
    odb::dbInst* dbInst = sink.getDbInst();
    if (dbInst == nullptr) return;
    int instTier = cts3dDb_->getInstTier(dbInst);
    if (instTier != tier) return;
    // FloatFix: (float) → (double) to prevent map key mismatch → CTS-0080
    sinkLocations.emplace_back((double) sink.getX() / wireSegmentUnit_,
                               (double) sink.getY() / wireSegmentUnit_);
    sinkInsts.emplace_back(&sink);
  });
}

void HTreeBuilder::initSecondLevelSinks(
    std::vector<std::pair<float, float>>& sinkLocations,
    std::vector<const ClockInst*>& sinkInsts)
{
  sinkLocations.clear();
  for (const auto& buf : topLevelSinksClustered_) {
    sinkLocations.emplace_back(buf.first, buf.second);
    const Point<double> bufPos(buf.first, buf.second);
    sinkInsts.emplace_back(mapLocationToSink_[bufPos]);
  }
}

void HTreeBuilder::computeBranchSinks(
    const LevelTopology& topology,
    const unsigned branchIdx,
    std::vector<std::pair<float, float>>& sinkLocations) const
{
  sinkLocations.clear();
  for (const Point<double>& point :
       topology.getBranchSinksLocations(branchIdx)) {
    sinkLocations.emplace_back(point.getX(), point.getY());
  }
}

void HTreeBuilder::refineBranchingPointsWithClustering(
    LevelTopology& topology,
    const unsigned level,
    const unsigned branchPtIdx1,
    const unsigned branchPtIdx2,
    const Point<double>& rootLocation,
    const std::vector<std::pair<float, float>>& sinks)
{
  CKMeans::Clustering clusteringEngine(
      sinks, rootLocation.getX(), rootLocation.getY(), logger_);

  // Pass per-sink LP targets to CKMeans.
  // When ckmeansTargetBeta_ > 0, CKMeans groups FFs with similar targets
  // into the same branch, enabling larger V31 wire shift differences.
  if (ckmeansTargetBeta_ > 0.0 && cts3dDb_ != nullptr
      && cts3dDb_->hasSkewTargets()) {
    std::vector<float> sinkTargets;
    sinkTargets.reserve(sinks.size());
    for (const auto& s : sinks) {
      const Point<double> sinkLoc(s.first, s.second);
      auto it = mapLocationToSink_.find(sinkLoc);
      if (it != mapLocationToSink_.end() && it->second != nullptr) {
        sinkTargets.push_back(
            static_cast<float>(getClockInstTarget(it->second)));
      } else {
        sinkTargets.push_back(0.0f);
      }
    }
    clusteringEngine.setSinkTargets(
        sinkTargets, static_cast<float>(ckmeansTargetBeta_));
  }

  Point<double>& branchPt1 = topology.getBranchingPoint(branchPtIdx1);
  Point<double>& branchPt2 = topology.getBranchingPoint(branchPtIdx2);

  std::vector<std::pair<float, float>> means;
  means.emplace_back(branchPt1.getX(), branchPt1.getY());
  means.emplace_back(branchPt2.getX(), branchPt2.getY());

  const unsigned cap
      = options_->getMaxFanout()
            ? (unsigned) (sinks.size() * 0.5)
            : (unsigned) (sinks.size() * options_->getClusteringCapacity());
  clusteringEngine.iterKmeans(
      1, means.size(), cap, 5, options_->getClusteringPower(), means);

  if (((int) options_->getNumStaticLayers() - (int) level) < 0) {
    branchPt1 = Point<double>(means[0].first, means[0].second);
    branchPt2 = Point<double>(means[1].first, means[1].second);
  }

  // Retrieve cluster assignments (used both for V31 wire shift and sink assignment)
  std::vector<std::vector<unsigned>> clusters;
  clusteringEngine.getClusters(clusters);

  // Useful-skew wire length adjustment (Fishburn 1990)
  // After k-means sets branch point positions, shift each branch point radially
  // from rootLocation based on its cluster's mean LP skew target.
  //
  // Physical meaning:
  //   shift_um = wireSkewScale_ * (T_cluster - T_global) / wireDelayPerUnit_
  //   T_cluster > T_global → branch farther from root → longer wire → later clock arrival
  //   T_cluster < T_global → branch closer to root  → shorter wire → earlier clock arrival
  //
  // This realizes the Fishburn LP-SAFETY schedule in the physical wire topology,
  // instead of patching with delay buffers after a balanced tree is built.
  if (wireSkewScale_ > 0.0 && cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()
      && clusters.size() >= 2) {
    double sumT[2] = {0.0, 0.0};
    int    cntT[2] = {0,   0  };
    for (int ci = 0; ci < 2; ++ci) {
      for (unsigned sinkIdx : clusters[ci]) {
        const Point<double> sinkLoc(sinks[sinkIdx].first, sinks[sinkIdx].second);
        auto it = mapLocationToSink_.find(sinkLoc);
        if (it != mapLocationToSink_.end() && it->second != nullptr) {
          // Unified target accessor (cluster buffer or individual FF)
          sumT[ci] += getClockInstTarget(it->second);
          cntT[ci]++;
        }
      }
    }
    const double meanT0     = (cntT[0] > 0) ? sumT[0] / cntT[0] : 0.0;
    const double meanT1     = (cntT[1] > 0) ? sumT[1] / cntT[1] : 0.0;
    const double globalMean = (meanT0 + meanT1) / 2.0;

    // Radially shift each branch point: positive shift = farther from root
    auto shiftBranchPoint = [&](Point<double>& bp, double shiftNs) {
      const double dx   = bp.getX() - rootLocation.getX();
      const double dy   = bp.getY() - rootLocation.getY();
      const double dist = std::sqrt(dx * dx + dy * dy);
      if (dist < 1e-6 || wireDelayPerUnit_ < 1e-12) {
        return;
      }
      const double shiftUm = wireSkewScale_ * shiftNs / wireDelayPerUnit_;
      // Widened clamp from 50% to 80% of branch length.
      // 50% limited wire-skew to ~10ps; LP requests 30-50ps differentials.
      const double maxShift    = 0.8 * dist;
      const double clampedShift = std::max(-maxShift, std::min(maxShift, shiftUm));
      const double scale = (dist + clampedShift) / dist;
      bp.setX(rootLocation.getX() + dx * scale);
      bp.setY(rootLocation.getY() + dy * scale);
    };

    shiftBranchPoint(branchPt1, meanT0 - globalMean);
    shiftBranchPoint(branchPt2, meanT1 - globalMean);

    // Log at debug level to avoid flooding; only log when shift is significant
    const double shift0Ps = wireSkewScale_ * (meanT0 - globalMean) / wireDelayPerUnit_ * 1000.0;
    const double shift1Ps = wireSkewScale_ * (meanT1 - globalMean) / wireDelayPerUnit_ * 1000.0;
    if (std::abs(shift0Ps) > 0.1 || std::abs(shift1Ps) > 0.1) {
      debugPrint(logger_, CTS, "skew_wire", 1,
                 "Level {} wire skew: cluster0 T={:+.1f}ps shift={:+.1f}um, "
                 "cluster1 T={:+.1f}ps shift={:+.1f}um",
                 level,
                 (meanT0 - globalMean) * 1000.0, shift0Ps / 1000.0 * wireDelayPerUnit_ / wireDelayPerUnit_,
                 (meanT1 - globalMean) * 1000.0, shift1Ps / 1000.0 * wireDelayPerUnit_ / wireDelayPerUnit_);
    }
  }
  // End wire-skew adjustment

  unsigned movedSinks = 0;
  const double errorFactor = 1.2;
  for (int clusterIdx = 0; clusterIdx < clusters.size(); ++clusterIdx) {
    for (int elementIdx = 0; elementIdx < clusters[clusterIdx].size();
         ++elementIdx) {
      const unsigned sinkIdx = clusters[clusterIdx][elementIdx];
      const Point<double> sinkLoc(sinks[sinkIdx].first, sinks[sinkIdx].second);
      const double dist = clusterIdx == 0 ? computeDist(branchPt1, sinkLoc)
                                          : computeDist(branchPt2, sinkLoc);
      const double distOther = clusterIdx == 0
                                   ? computeDist(branchPt2, sinkLoc)
                                   : computeDist(branchPt1, sinkLoc);
      if (clusterIdx == 0) {
        topology.addSinkToBranch(branchPtIdx1, sinkLoc);
      } else {
        topology.addSinkToBranch(branchPtIdx2, sinkLoc);
      }

      if (dist >= distOther * errorFactor) {
        movedSinks++;
      }
    }
  }

  if (movedSinks > 0) {
    debugPrint(logger_,
               CTS,
               "clustering",
               1,
               " Out of {} sinks, {} sinks closer to other cluster.",
               sinks.size(),
               movedSinks);
  }
}

// Unified LP skew target accessor for any ClockInst.
// After preSinkClustering(), mapLocationToSink_ entries for multi-FF clusters
// point to cluster buffers (clkbuf_leaf_X), NOT individual FFs.
// Cluster buffer names are not in the LP CSV, so direct name lookup returns 0.
// This helper checks clusterBufTarget_ first (set at clustering time from FF members),
// then falls back to stripping "/clk_pin" suffix for individual FF lookup.
double HTreeBuilder::getClockInstTarget(ClockInst* inst) const
{
  if (inst == nullptr || cts3dDb_ == nullptr || !cts3dDb_->hasSkewTargets()) {
    return 0.0;
  }
  // Cluster buffer path: use precomputed mean of member FF targets
  auto cit = clusterBufTarget_.find(inst);
  if (cit != clusterBufTarget_.end()) {
    return cit->second;
  }
  // Individual FF path: ClockInst name = "inst_name/clk_pin" → strip suffix
  std::string nm = inst->getName();
  const auto sp = nm.rfind('/');
  if (sp != std::string::npos) nm = nm.substr(0, sp);
  return cts3dDb_->getSkewTarget(nm);
}

// Step 4: Compute mean arrival target per branch.
// Extended to ALL levels (not just leaf).
// This determines which branches need extra delay buffers at each level.
void HTreeBuilder::computeBranchDelayTargets()
{
  branchDelayTargets_.clear();
  globalMeanTarget_ = 0.0;
  perLevelGlobalMean_.clear();

  if (topologyForEachLevel_.empty() || cts3dDb_ == nullptr
      || !cts3dDb_->hasSkewTargets()) {
    return;
  }

  // Debug: per-branch LP target statistics
  struct BranchDetail {
    int level;
    unsigned idx;
    int count;
    double mean, bmin, bmax, stddev;
  };
  std::vector<BranchDetail> branchDetails;

  // Compute per-branch mean target at ALL levels (not just leaf).
  // At each level, a branch's mean target = average of all sink targets
  // reachable from that branch point downward.
  for (int levelIdx = 0;
       levelIdx < static_cast<int>(topologyForEachLevel_.size());
       ++levelIdx) {
    LevelTopology& topology = topologyForEachLevel_[levelIdx];
    double levelTargetSum = 0.0;
    int levelSinkCount = 0;

    topology.forEachBranchingPoint([&](unsigned idx, Point<double>) {
      const auto& sinkLocs = topology.getBranchSinksLocations(idx);
      if (sinkLocs.empty()) {
        return;
      }

      double branchSum = 0.0, branchSumSq = 0.0;
      double branchMin = std::numeric_limits<double>::max();
      double branchMax = -std::numeric_limits<double>::max();
      int branchCount = 0;

      for (const auto& loc : sinkLocs) {
        // V58 Fix B: epsilon-tolerant lookup (was exact find → skipped
        // singleton macro clusters whose float→double key didn't match)
        ClockInst* sink = findSinkEps(loc);
        if (sink == nullptr) {
          continue;
        }
        const double tgt = getClockInstTarget(sink);
        branchSum   += tgt;
        branchSumSq += tgt * tgt;
        branchMin    = std::min(branchMin, tgt);
        branchMax    = std::max(branchMax, tgt);
        ++branchCount;
      }

      if (branchCount > 0) {
        const double meanTarget = branchSum / branchCount;
        const double var = (branchSumSq / branchCount) - meanTarget * meanTarget;
        const double stddev = (var > 0.0) ? std::sqrt(var) : 0.0;
        branchDelayTargets_[{levelIdx, idx}] = meanTarget;
        levelTargetSum += branchSum;
        levelSinkCount += branchCount;
        branchDetails.push_back(
            {levelIdx, idx, branchCount, meanTarget,
             branchMin, branchMax, stddev});
      }
    });

    if (levelSinkCount > 0) {
      perLevelGlobalMean_[levelIdx] = levelTargetSum / levelSinkCount;
    }
  }

  // Global mean across all sinks (backward compatible with V32c leaf-only)
  const int leafLevelIdx = static_cast<int>(topologyForEachLevel_.size()) - 1;
  auto it = perLevelGlobalMean_.find(leafLevelIdx);
  globalMeanTarget_ = (it != perLevelGlobalMean_.end()) ? it->second : 0.0;

  // Log per-level statistics
  for (int levelIdx = 0;
       levelIdx < static_cast<int>(topologyForEachLevel_.size());
       ++levelIdx) {
    int posCount = 0, negCount = 0;
    int branchCount = 0;
    const double levelMean
        = perLevelGlobalMean_.count(levelIdx)
            ? perLevelGlobalMean_[levelIdx] : 0.0;
    for (const auto& [key, meanTarget] : branchDelayTargets_) {
      if (key.first != levelIdx) continue;
      ++branchCount;
      const double delta = meanTarget - levelMean;
      if (delta > delayTargetThreshold_) ++posCount;
      else if (delta < -delayTargetThreshold_) ++negCount;
    }
    logger_->info(CTS, 401,
                  "Per-branch delay targets L{}: {} branches, "
                  "levelMean={:.4f}ns, {} need extra delay, {} need less delay",
                  levelIdx + 1, branchCount, levelMean, posCount, negCount);
  }

  // Per-branch LP target distribution detail (leaf level only)
  for (const auto& d : branchDetails) {
    if (d.level != leafLevelIdx) continue;
    // globalMean fix: absolute LP target for planned_N diagnostic
    const double delta = d.mean;
    int plannedN = 0;
    if (delta > 0.0 && singleBufDelay_ > 1e-9) {
      plannedN = std::max(0, std::min(
          static_cast<int>(std::round(delta / singleBufDelay_)),
          maxLeafDelayBufs_));
    }
    debugPrint(logger_, CTS, "HTree", 1,
                  "  Leaf branch {:2d}: sinks={}, "
                  "LP=[min={:.4f} mean={:.4f} max={:.4f} spread={:.4f} "
                  "stddev={:.4f}]ns, delta={:.4f}ns -> planned_N={}",
                  d.idx, d.count,
                  d.bmin, d.mean, d.bmax, d.bmax - d.bmin,
                  d.stddev, delta, plannedN);
  }
}

// Load per-FF hold budgets from timing graph CSV.
// For each capture FF (to_ff column), compute min(slack_min_ns) across all
// edges. This is the FF's hold budget: max additional delay we can add
// to its clock arrival without violating hold.
void HTreeBuilder::loadPerFfHoldBudgets(const std::string& timingGraphCsv)
{
  perFfHoldBudget_.clear();
  std::ifstream file(timingGraphCsv);
  if (!file.is_open()) {
    logger_->warn(CTS, 422,
                  "Cannot open timing graph for hold budgets: {}",
                  timingGraphCsv);
    return;
  }

  std::string line;
  std::getline(file, line);  // skip header

  int edgeCount = 0;
  while (std::getline(file, line)) {
    std::istringstream ss(line);
    std::string from_ff, to_ff, slack_max_str, slack_min_str;
    // CSV format: from_ff,to_ff,slack_max_ns,slack_min_ns,...
    std::getline(ss, from_ff, ',');
    std::getline(ss, to_ff, ',');
    std::getline(ss, slack_max_str, ',');
    std::getline(ss, slack_min_str, ',');

    if (to_ff.empty() || slack_min_str.empty()) continue;

    try {
      const double holdSlack = std::stod(slack_min_str);
      auto it = perFfHoldBudget_.find(to_ff);
      if (it == perFfHoldBudget_.end() || holdSlack < it->second) {
        perFfHoldBudget_[to_ff] = holdSlack;
      }
      ++edgeCount;
    } catch (...) {
      // skip malformed lines
    }
  }

  logger_->info(CTS, 423,
                "Hold budgets loaded: {} capture FFs from {} edges",
                perFfHoldBudget_.size(), edgeCount);
}

// Load IO edge hold budgets from IO timing CSV.
// For PI_TO_FF rows, the hold slack represents the budget for adding delay
// to the capture FF's clock arrival. Merged into perFfHoldBudget_ (take min).
void HTreeBuilder::loadIoHoldBudgets(const std::string& ioCsv)
{
  std::ifstream file(ioCsv);
  if (!file.is_open()) {
    logger_->warn(CTS, 429,
                  "Cannot open IO timing CSV for hold budgets: {}", ioCsv);
    return;
  }

  std::string line;
  std::getline(file, line);  // skip header: edge_type,port_name,ff_name,...

  int ioEdgeCount = 0;
  int mergedCount = 0;
  int newCount = 0;
  while (std::getline(file, line)) {
    std::istringstream ss(line);
    std::string edgeType, portName, ffName, slackSetup, slackHold;
    // CSV: edge_type,port_name,ff_name,slack_setup_ns,slack_hold_ns
    std::getline(ss, edgeType, ',');
    std::getline(ss, portName, ',');
    std::getline(ss, ffName, ',');
    std::getline(ss, slackSetup, ',');
    std::getline(ss, slackHold, ',');

    // Only PI_TO_FF hold matters: adding relay to capture FF hurts hold.
    // FF_TO_PO: adding relay to launch FF helps hold (not a concern).
    if (edgeType != "PI_TO_FF") continue;
    if (ffName.empty() || slackHold.empty()) continue;

    try {
      const double holdSlack = std::stod(slackHold);
      auto it = perFfHoldBudget_.find(ffName);
      if (it == perFfHoldBudget_.end()) {
        perFfHoldBudget_[ffName] = holdSlack;
        ++newCount;
      } else if (holdSlack < it->second) {
        it->second = holdSlack;
        ++mergedCount;
      }
      ++ioEdgeCount;
    } catch (...) {
      // skip malformed lines
    }
  }

  logger_->info(CTS, 430,
                "IO hold budgets: {} PI→FF edges, {} new FFs, "
                "{} tightened (total hold budget FFs: {})",
                ioEdgeCount, newCount, mergedCount,
                perFfHoldBudget_.size());
}

void HTreeBuilder::createClockSubNets()
{
  Point<double> legalCenter{0.0, 0.0};
  ClockInst* rootBufferPtr;
  ClockSubNet* rootSubNetPtr;

  if (sharedRootBuffer_ != nullptr) {
    // Per-tier mode: reuse the shared root buffer created by createSharedRoot()
    rootBufferPtr = sharedRootBuffer_;
    rootSubNetPtr = sharedRootSubNet_;
    legalCenter = sharedRootLocation_;
    logger_->info(utl::CTS, 857,
                  "Per-tier sub-tree (tier prefix='{}') using shared root "
                  "at ({:.2f}, {:.2f})",
                  tierPrefix_, legalCenter.getX(), legalCenter.getY());
  } else {
    // Unified mode: create root buffer here (original behavior)
    Point<double> center = sinkRegion_.getCenter();
    // V58_CL2: CTS_FORCE_ROOT_TIER overrides dominant-tier selection (same as createSharedRoot)
    int root_tier = cts3dDb_->getDominantTierFromClock(clock_);
    if (const char* e = std::getenv("CTS_FORCE_ROOT_TIER")) {
      int forced = std::atoi(e);
      if (forced >= 0 && forced != root_tier) {
        logger_->info(CTS, 701,
            "V58_CL2: Root tier forced {} -> {} (CTS_FORCE_ROOT_TIER, unified)",
            root_tier, forced);
        root_tier = forced;
      }
    }
    const std::string root_buffer
        = cts3dDb_->getBufferForTier(options_->getRootBuffer(), root_tier);
    logger_->info(utl::CTS, 315,
                  "3D-CTS root buffer: tier={}, buffer={}",
                  root_tier, root_buffer);
    legalCenter = legalizeOneBuffer(center, root_buffer);
    sinkRegion_.setCenter(legalCenter);
    commitMoveLoc(center, legalCenter);
    const int centerX = legalCenter.getX() * wireSegmentUnit_;
    const int centerY = legalCenter.getY() * wireSegmentUnit_;

    const std::string rootBufName = tierPrefix_ + "clkbuf_0";
    const std::string rootNetName = tierPrefix_ + "clknet_0";

    ClockInst& rootBuffer
        = clock_.addClockBuffer(rootBufName, root_buffer, centerX, centerY);

    if (topBufferName_.empty()) {
      topBufferName_ = rootBuffer.getName();
    }

    // clang-format off
    if (center != legalCenter) {
      debugPrint(logger_, CTS, "legalizer", 2, "createClockSubNets: "
                 "root {}: {} => {}", rootBufName, center, legalCenter);
    } else {
      debugPrint(logger_, CTS, "legalizer", 2, "createClockSubNets: "
                 "root {}: {}", rootBufName, center);
    }
    // clang-format on

    addTreeLevelBuffer(&rootBuffer);
    ClockSubNet& rootClockSubNet = clock_.addSubNet(rootNetName);
    rootClockSubNet.setTreeLevel(0);
    rootClockSubNet.addInst(rootBuffer);
    treeBufLevels_++;

    rootBufferPtr = &rootBuffer;
    rootSubNetPtr = &rootClockSubNet;
  }

  // First level...
  LevelTopology& topLevelTopology = topologyForEachLevel_[0];
  bool isFirstPoint = true;
  topLevelTopology.forEachBranchingPoint([&](unsigned idx,
                                             Point<double> branchPoint) {
    // If the branch point has no sinks that will be connected to
    // it don't create a clock sub net for it
    if (topLevelTopology.getBranchSinksLocations(idx).empty()) {
      return;
    }
    // Replaced with Cts3DDatabase calls
    const int branch_tier = cts3dDb_->getDominantTier(
        topLevelTopology.getBranchSinksLocations(idx), mapLocationToSink_);
    const std::string branch_root_buffer = cts3dDb_->getBufferForTier(
        options_->getRootBuffer(), branch_tier);
    logger_->info(utl::CTS, 316,
                  "3D-CTS branch L1 idx={}: tier={}, buffer={}",
                  idx, branch_tier, branch_root_buffer);

    // Cross-tier HB delay consideration
    // Elmore 50% delay: t_via = 0.693 * R_HB * C_HB (negligible: ~0 with C=0)
    // In shared-root mode, root_tier comes from the shared root buffer's tier.
    const int root_tier = cts3dDb_->getDominantTierFromClock(clock_);
    if (branch_tier != root_tier) {
      const double hb_res = cts3dDb_->getHbtResistance();
      const double hb_cap = cts3dDb_->getHbtCapacitance();
      const double hb_delay = 0.693 * hb_res * hb_cap;

      const double wire_res_per_unit = cts3dDb_->getResPerDBU(branch_tier) * wireSegmentUnit_;
      const double wire_cap_per_unit = cts3dDb_->getCapPerDBU(branch_tier) * wireSegmentUnit_;
      const double hb_equivalent_dist = cts3dDb_->getHbtEquivalentDistance(
          wire_res_per_unit, wire_cap_per_unit);

      logger_->warn(CTS, 362,
                    "Cross-tier branch L1: branch_tier={}, root_tier={}, "
                    "HB_delay={:.3e}s, equiv_dist={:.1f}um",
                    branch_tier, root_tier, hb_delay, hb_equivalent_dist);
    }

    Point<double> legalBranchPoint
        = legalizeOneBuffer(branchPoint, branch_root_buffer);
    commitMoveLoc(branchPoint, legalBranchPoint);

    // clang-format off
    if (branchPoint != legalBranchPoint) {
      debugPrint(logger_, CTS, "legalizer", 2,
		 "createClockSubNets level 1 clk_buf_1_{}_ : {} => {}",
		 std::to_string(idx), branchPoint, legalBranchPoint);
    } else {
      debugPrint(logger_, CTS, "legalizer", 2, 
		 "createClockSubNets level 1 clk_buf_1_{}_ : {}",
		 std::to_string(idx), branchPoint);
    }
    // clang-format on

    SegmentBuilder builder(tierPrefix_ + "clkbuf_1_" + std::to_string(idx) + "_",
                           tierPrefix_ + "clknet_1_" + std::to_string(idx) + "_",
                           legalCenter,  // center may have moved, don't use
                                         // sinkRegion_.getCenter()
                           legalBranchPoint,
                           topLevelTopology.getWireSegments(),
                           clock_,
                           *rootSubNetPtr,
                           *techChar_,
                           wireSegmentUnit_,
                           this,
                           db_,
                           branch_tier);
    if (!options_->getTreeBuffer().empty()) {
      // Replaced mapBufferMasterToTier with Cts3DDatabase
      const std::string tree_buffer = cts3dDb_->getBufferForTier(
          options_->getTreeBuffer(), branch_tier);
      builder.build(tree_buffer);
    } else {
      builder.build();
    }
    if (topologyForEachLevel_.size() == 1) {
      builder.forceBufferInSegment(branch_root_buffer);
    }
    if (isFirstPoint) {
      treeBufLevels_ += builder.getNumBufferLevels();
      isFirstPoint = false;
    }

    // Insert intermediate delay buffers at Level 1.
    // If this branch's subtree mean target > level globalMean, add N delay
    // buffers in series to shift arrival for the entire subtree.
    ClockSubNet* branchDrivingSub = builder.getDrivingSubNet();
    // gate with enableMidDelayBufs_ (default off in V39)
    if (enableMidDelayBufs_ && cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()
        && singleBufDelay_ > 1e-9) {
      auto tgtIt = branchDelayTargets_.find({0, idx});  // level 0
      const double levelMean
          = perLevelGlobalMean_.count(0) ? perLevelGlobalMean_[0] : 0.0;
      if (tgtIt != branchDelayTargets_.end()) {
        const double delta = tgtIt->second - levelMean;
        int midN = 0;
        if (delta > 0.0) {
          midN = static_cast<int>(std::round(delta / singleBufDelay_));
          midN = std::max(0, std::min(midN, maxLeafDelayBufs_));
        }
        if (midN > 0) {
          const int bx = legalBranchPoint.getX() * wireSegmentUnit_;
          const int by = legalBranchPoint.getY() * wireSegmentUnit_;
          const std::string delayBufMaster = branch_root_buffer;
          ClockSubNet* curSub = branchDrivingSub;
          for (int b = 0; b < midN; ++b) {
            const std::string bufName = tierPrefix_ + "clkbuf_mid_dly_1_"
                + std::to_string(idx) + "_" + std::to_string(b);
            ClockInst& buf = clock_.addClockBuffer(
                bufName, delayBufMaster, bx, by);
            cts3dDb_->setClockInstTier(buf, branch_tier);
            addTreeLevelBuffer(&buf);
            curSub->addInst(buf);
            ClockSubNet& nextSub = clock_.addSubNet(
                tierPrefix_ + "clknet_mid_dly_1_" + std::to_string(idx)
                + "_" + std::to_string(b));
            nextSub.addInst(buf);
            curSub = &nextSub;
          }
          branchDrivingSub = curSub;
          logger_->info(CTS, 411,
                        "Intermediate delay L1 branch {}: "
                        "target={:.4f}ns, levelMean={:.4f}ns, "
                        "delta={:.4f}ns, N={}",
                        idx, tgtIt->second, levelMean, delta, midN);
        }
      }
    }
    topLevelTopology.setBranchDrivingSubNet(idx, *branchDrivingSub);
  });

  // Others...
  for (int levelIdx = 1; levelIdx < topologyForEachLevel_.size(); ++levelIdx) {
    LevelTopology& topology = topologyForEachLevel_[levelIdx];
    isFirstPoint = true;
    topology.forEachBranchingPoint([&](unsigned idx,
                                       Point<double> branchPoint) {
      // If the branch point has no sinks that will be connected
      // to it don't create a clock sub net for it
      if (topology.getBranchSinksLocations(idx).empty()) {
        return;
      }
      unsigned parentIdx = topology.getBranchingPointParentIdx(idx);
      LevelTopology& parentTopology = topologyForEachLevel_[levelIdx - 1];
      Point<double> parentPoint = parentTopology.getBranchingPoint(parentIdx);

      // Replaced with Cts3DDatabase calls
      const int branch_tier = cts3dDb_->getDominantTier(
          topology.getBranchSinksLocations(idx), mapLocationToSink_);
      const std::string branch_root_buffer = cts3dDb_->getBufferForTier(
          options_->getRootBuffer(), branch_tier);
      debugPrint(logger_, CTS, "HTree", 1,
                    "3D-CTS branch L{} idx={}: tier={}, buffer={}",
                    levelIdx+1, idx, branch_tier, branch_root_buffer);

      // Cross-tier HB delay consideration
      // Elmore 50% delay: t_via = 0.693 * R_HB * C_HB (negligible: ~0 with C=0)
      const int parent_tier = cts3dDb_->getDominantTier(
          parentTopology.getBranchSinksLocations(parentIdx), mapLocationToSink_);
      if (branch_tier != parent_tier) {
        const double hb_res = cts3dDb_->getHbtResistance();
        const double hb_cap = cts3dDb_->getHbtCapacitance();
        const double hb_delay = 0.693 * hb_res * hb_cap;

        const double wire_res_per_unit = cts3dDb_->getResPerDBU(branch_tier) * wireSegmentUnit_;
        const double wire_cap_per_unit = cts3dDb_->getCapPerDBU(branch_tier) * wireSegmentUnit_;
        const double hb_equivalent_dist = cts3dDb_->getHbtEquivalentDistance(
            wire_res_per_unit, wire_cap_per_unit);

        logger_->warn(CTS, 363,
                      "Cross-tier branch L{}: branch_tier={}, parent_tier={}, "
                      "HB_delay={:.3e}s, equiv_dist={:.1f}um",
                      levelIdx+1, branch_tier, parent_tier, hb_delay, hb_equivalent_dist);
      }

      Point<double> legalBranchPoint
          = legalizeOneBuffer(branchPoint, branch_root_buffer);
      commitMoveLoc(branchPoint, legalBranchPoint);

      // clang-format off
      if (branchPoint != legalBranchPoint) {
	debugPrint(logger_, CTS, "legalizer", 2, "createClockSubNets level {} "
		   "{} : {} => {}", levelIdx+1, tierPrefix_ + "clkbuf_" +
		   std::to_string(levelIdx+1) + "_" + std::to_string(idx) + "_",
		   branchPoint, legalBranchPoint);
      } else {
	debugPrint(logger_, CTS, "legalizer", 2, "createClockSubNets level {} "
		   "{} : {}", levelIdx+1, tierPrefix_ + "clkbuf_" +
		   std::to_string(levelIdx+1) + "_" + std::to_string(idx) + "_",
		   branchPoint);
      }
      // clang-format on

      SegmentBuilder builder(tierPrefix_ + "clkbuf_" + std::to_string(levelIdx + 1) + "_"
                                 + std::to_string(idx) + "_",
                             tierPrefix_ + "clknet_" + std::to_string(levelIdx + 1) + "_"
                                 + std::to_string(idx) + "_",
                             parentPoint,
                             legalBranchPoint,
                             topology.getWireSegments(),
                             clock_,
                             *parentTopology.getBranchDrivingSubNet(parentIdx),
                             *techChar_,
                             wireSegmentUnit_,
                             this,
                             db_,
                             branch_tier);

      // Set clock tree level the first time only.
      if (builder.getDrivingSubNet()->getTreeLevel() < 0) {
        builder.getDrivingSubNet()->setTreeLevel(levelIdx);
      }

      if (!options_->getTreeBuffer().empty()) {
        // Replaced mapBufferMasterToTier with Cts3DDatabase
        const std::string tree_buffer = cts3dDb_->getBufferForTier(
            options_->getTreeBuffer(), branch_tier);
        builder.build(tree_buffer);
      } else {
        builder.build();
      }
      if (levelIdx == topologyForEachLevel_.size() - 1) {
        builder.forceBufferInSegment(branch_root_buffer);
      }
      if (isFirstPoint) {
        treeBufLevels_ += builder.getNumBufferLevels();
        isFirstPoint = false;
      }

      // Insert intermediate delay buffers at Level 2+.
      // Skip the leaf level (handled by existing N-chain code below).
      const int leafLevel
          = static_cast<int>(topologyForEachLevel_.size()) - 1;
      ClockSubNet* branchDrivingSub = builder.getDrivingSubNet();
      // gate with enableMidDelayBufs_ (default off in V39)
      if (enableMidDelayBufs_ && levelIdx < leafLevel
          && cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()
          && singleBufDelay_ > 1e-9) {
        auto tgtIt = branchDelayTargets_.find({levelIdx, idx});
        const double levelMean
            = perLevelGlobalMean_.count(levelIdx)
                ? perLevelGlobalMean_[levelIdx] : 0.0;
        if (tgtIt != branchDelayTargets_.end()) {
          const double delta = tgtIt->second - levelMean;
          int midN = 0;
          if (delta > 0.0) {
            midN = static_cast<int>(std::round(delta / singleBufDelay_));
            midN = std::max(0, std::min(midN, maxLeafDelayBufs_));
          }
          if (midN > 0) {
            const int bx = legalBranchPoint.getX() * wireSegmentUnit_;
            const int by = legalBranchPoint.getY() * wireSegmentUnit_;
            const std::string delayBufMaster = branch_root_buffer;
            ClockSubNet* curSub = branchDrivingSub;
            for (int b = 0; b < midN; ++b) {
              const std::string bufName = tierPrefix_ + "clkbuf_mid_dly_"
                  + std::to_string(levelIdx + 1) + "_"
                  + std::to_string(idx) + "_" + std::to_string(b);
              ClockInst& buf = clock_.addClockBuffer(
                  bufName, delayBufMaster, bx, by);
              cts3dDb_->setClockInstTier(buf, branch_tier);
              addTreeLevelBuffer(&buf);
              curSub->addInst(buf);
              ClockSubNet& nextSub = clock_.addSubNet(
                  tierPrefix_ + "clknet_mid_dly_" + std::to_string(levelIdx + 1)
                  + "_" + std::to_string(idx) + "_" + std::to_string(b));
              nextSub.addInst(buf);
              curSub = &nextSub;
            }
            branchDrivingSub = curSub;
            logger_->info(CTS, 412,
                          "Intermediate delay L{} branch {}: "
                          "target={:.4f}ns, levelMean={:.4f}ns, "
                          "delta={:.4f}ns, N={}",
                          levelIdx + 1, idx, tgtIt->second,
                          levelMean, delta, midN);
          }
        }
      }
      topology.setBranchDrivingSubNet(idx, *branchDrivingSub);
    });
  }

  LevelTopology& leafTopology = topologyForEachLevel_.back();
  // leafLevelIdx removed: no longer needed (branchDelayTargets_ lookup replaced by inline absTgt)
  unsigned numSinks = 0;
  unsigned numDelayBufs = 0;

  // Debug: LP-CTS gap tracking accumulators (reported in CTS-408 after loop).
  // totalLpDelay   = sum(branchMean * sinkCount) over all leaf branches
  // totalDelivered = sum(N * singleBufDelay * sinkCount) for fired branches only
  // Coverage = totalDelivered / totalLpDelay * 100%: how much LP plan was implemented.
  double totalLpDelay_debug   = 0.0;
  double totalDelivered_debug = 0.0;
  int totalBranchSinks_debug  = 0;
  int totalBranches_debug     = 0;
  int firedBranches_debug     = 0;

  leafTopology.forEachBranchingPoint(
      [&](unsigned idx, Point<double> branchPoint) {
        ClockSubNet* subNet = leafTopology.getBranchDrivingSubNet(idx);
        // If no clock sub net was created for a leaf branch point no sinks
        // connect to it, so just skip.
        if (subNet == nullptr) {
          return;
        }

        // Step 4: Insert delay buffers if branch LP target > globalMean
        // N-buffer chain based on ABSOLUTE LP-TNS target.
        //   Problem (found after V32b run): Absolute targeting causes hold violations.
        //   Branch differentials up to 36-48ps (N=1 vs N=4) violated hold on short paths.
        //   LP guaranteed hold only for per-FF delays; branch-mean averaging breaks this.
        // DELTA-BASED N-buffer chain (relative to globalMean).
        //   N = round(max(0, branchMean - globalMean) / singleBufDelay_)
        //   Limits max inter-branch differential to (maxBranchDelta - 0) ≈ 12-16ps (1-2 bufs).
        //   LP-TNS result (asap7/aes): globalMean=26.6ps, max delta=15.8ps → N≤1 per branch.
        //   Chain: subNet → buf_0 → sub_0 → ... → FFs
        // Debug: branchCnt hoisted to lambda scope for CTS-408 gap tracking
        double absTgt = 0.0;
        int branchCnt = 0;
        if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()
            && singleBufDelay_ > 1e-9) {
          const auto& bSinkLocs = leafTopology.getBranchSinksLocations(idx);
          double branchSum = 0.0;
          for (const auto& bLoc : bSinkLocs) {
            // V58 Fix B: epsilon-tolerant lookup for branch target calculation
            ClockInst* bSink = findSinkEps(bLoc);
            if (bSink == nullptr) {
              continue;
            }
            // Unified target accessor (cluster buffer or individual FF)
            branchSum += getClockInstTarget(bSink);
            ++branchCnt;
          }
          if (branchCnt > 0) {
            absTgt = branchSum / branchCnt;
          }
        }

        // Grouped Delay Chain.
        // Instead of per-FF relay (1 chain per FF → ~192 nets → GRT-0183),
        // group FFs within each leaf cluster by quantized LP target delta.
        // Group 0 (delta < min): direct to leaf buffer.
        // Group k>0: shared chain of k delay buffers → all group FFs.
        // ~24 group nets (8 branches × ~3 groups) → no GRT-0183.
        //
        // Structure: branch → leaf_buffer → FFs (group 0, direct)
        //                                 → [grp_buf_0] → FFs (group 1)
        //                                 → [grp_buf_0] → [grp_buf_1] → FFs (group 2)
        if (enableGroupedDelay_ && cts3dDb_ != nullptr
            && cts3dDb_->hasSkewTargets() && perFfBufDelay_ > 1e-9) {
          // static counter for globally unique buffer names.
          // Multiple sub-trees (register tree, macro tree, etc.) each call
          // this code with idx=0,1,... → name collision without global ID.
          static int grpGlobalId = 0;
          const int myGrpId = grpGlobalId++;

          subNet->setLeafLevel(true);
          int branchGroupedFFs = 0, branchGroupBufs = 0;
          int branchDirectFFs = 0, branchHoldSkip = 0;

          const int branch_tier = cts3dDb_->getDominantTier(
              leafTopology.getBranchSinksLocations(idx), mapLocationToSink_);
          const std::string bufMaster = cts3dDb_->getBufferForTier(
              options_->getRootBuffer(), branch_tier);

          // Helper: get FF name (strip /CLK suffix)
          auto getFFName = [](ClockInst* ff) -> std::string {
            std::string nm = ff->getName();
            const auto sp = nm.rfind('/');
            if (sp != std::string::npos) nm = nm.substr(0, sp);
            return nm;
          };

          // Helper: get per-FF target delta and hold-safe depth
          struct FfGroupInfo {
            ClockInst* inst;
            double delta;
            int depth;  // quantized, hold-capped
          };

          // Collect all FFs in this branch with their group assignment.
          // When enableClusterUniform_=true, all FFs in the cluster share
          // a single depth k_c computed from the cluster median LP target delta.
          // This creates exactly 1 chain per cluster (vs up to MAX_DEPTH chains
          // in V49), reducing buffer count by ~60% with comparable timing.
          // Hold safety: k_c is capped by the minimum hold budget across all FFs.
          auto collectFFs = [&](ClockInst* driverInst, ClockSubNet* driverSubNet,
                                std::vector<FfGroupInfo>& ffInfos,
                                std::vector<ClockInst*>& memberFFs) {
            if (enableClusterUniform_) {
              // Cluster-Uniform mode — compute median delta, assign to all FFs.
              // Pass 1: collect deltas and hold budgets
              std::vector<double> deltas;
              deltas.reserve(memberFFs.size());
              for (ClockInst* ff : memberFFs) {
                std::string ffName = getFFName(ff);
                // globalMean fix: absolute LP target for TAP delivery
                double ffDelta = cts3dDb_->getSkewTarget(ffName);
                // Debug: trace per-FF target lookup (cluster-uniform path)
                if (ffDelta > 0.01) {
                  logger_->info(CTS, 711,
                      "collectFFs cluster-uniform: '{}' target={:.4f}ns",
                      ffName, ffDelta);
                }
                deltas.push_back(ffDelta);
              }
              // Compute median delta for the cluster
              std::vector<double> sorted = deltas;
              std::sort(sorted.begin(), sorted.end());
              double medianDelta = sorted[sorted.size() / 2];

              // Quantize cluster depth from median
              int clusterDepth = 0;
              if (medianDelta > groupedDelayMinDelta_) {
                clusterDepth = static_cast<int>(
                    std::round(medianDelta / perFfBufDelay_));
                clusterDepth = std::max(
                    0, std::min(clusterDepth, groupedDelayMaxDepth_));
              }

              // Removed redundant C++ holdSkip.
              // LP already enforces hold via hold edges + PI-hold-clip + sigma guard.
              // C++ holdSkip was double-guarding with a more conservative local metric,
              // blocking 60-70% of grouped delay delivery. See CLAUDE.md for details.

              // Pass 2: assign uniform depth to all FFs
              for (size_t i = 0; i < memberFFs.size(); ++i) {
                // V58 Fix C: Macro sinks (SRAM) get depth=0 — tight hold margin.
                int ffDepth = clusterDepth;
                odb::dbInst* dInst = memberFFs[i]->getDbInst();
                if (isMacroBlockInst(dInst)) {
                  ffDepth = 0;
                }
                ffInfos.push_back({memberFFs[i], deltas[i], ffDepth});
              }
            } else {
              // Per-FF depth (original behavior)
              for (ClockInst* ff : memberFFs) {
                std::string ffName = getFFName(ff);
                double ffTarget = cts3dDb_->getSkewTarget(ffName);
                // Debug: trace per-FF target lookup to diagnose delivery gap
                // (3058 non-zero LP targets but only ~29 grouped FFs)
                if (ffTarget > 0.01) {
                  debugPrint(logger_, CTS, "HTree", 1,
                      "collectFFs per-FF: '{}' target={:.4f}ns depth_will_be={}",
                      ffName, ffTarget,
                      static_cast<int>(std::round(ffTarget / perFfBufDelay_)));
                }
                // globalMean fix: absolute LP target
                double ffDelta = ffTarget;

                // Quantize to buffer delay steps
                int depth = 0;
                if (ffDelta > groupedDelayMinDelta_) {
                  depth = static_cast<int>(
                      std::round(ffDelta / perFfBufDelay_));
                  depth = std::max(0, std::min(depth, groupedDelayMaxDepth_));
                }

                // Removed redundant holdSkip (same as above).
                // V58 Fix C: Macro sinks (SRAM) get depth=0.
                odb::dbInst* dInst = ff->getDbInst();
                if (isMacroBlockInst(dInst)) {
                  depth = 0;
                }
                ffInfos.push_back({ff, ffDelta, depth});
              }
            }
          };

          // V58_CL: Hold-budget-aware TAP depth cap.
          // After collectFFs assigns depth from LP targets, check each FF's
          // hold budget from perFfHoldBudget_ (loaded from timing graph).
          // If depth * perFfBufDelay_ would exceed the FF's hold budget,
          // cap depth to floor(holdBudget / perFfBufDelay_).
          // This protects FFs where LP missed hold constraints (pruned edges)
          // but the actual timing graph has tight PI→FF hold paths.
          // CTS_ENABLE_HOLD_BUDGET_CAP=1 (default ON when hold budgets loaded).
          bool enableHoldBudgetCap = !perFfHoldBudget_.empty();
          if (const char* e = std::getenv("CTS_ENABLE_HOLD_BUDGET_CAP")) {
            enableHoldBudgetCap = (std::atoi(e) != 0);
          }
          auto applyHoldBudgetCap = [&](std::vector<FfGroupInfo>& ffInfos) {
            if (!enableHoldBudgetCap || perFfHoldBudget_.empty()) return;
            int nCapped = 0;
            for (auto& fi : ffInfos) {
              if (fi.depth <= 0) continue;
              std::string ffName = getFFName(fi.inst);
              auto it = perFfHoldBudget_.find(ffName);
              if (it == perFfHoldBudget_.end()) continue;
              double holdBudgetNs = it->second;
              // Max depth that fits within hold budget
              int maxSafeDepth = static_cast<int>(
                  std::floor(holdBudgetNs / perFfBufDelay_));
              if (maxSafeDepth < 0) maxSafeDepth = 0;
              if (fi.depth > maxSafeDepth) {
                fi.depth = maxSafeDepth;
                ++nCapped;
              }
            }
            if (nCapped > 0) {
              logger_->info(CTS, 698,
                  "V58_CL: Hold budget cap reduced TAP depth for {} FFs "
                  "(bufDelay={:.3f}ns)",
                  nCapped, perFfBufDelay_);
            }
          };

          // Build grouped delay chains for collected FFs
          // No local numSinks here. Uses outer numSinks (line 3002)
          // so buildGroupChains [&] captures the correct per-sink counter.
          // Previously an inner numSinks was declared AFTER this lambda,
          // so the lambda captured the outer numSinks (always 0) → duplicate
          // buffer names across cluster buffers in same branch → CTS-0499.
          auto buildGroupChains = [&](ClockInst* driverInst,
                                      ClockSubNet* driverSubNet,
                                      std::vector<FfGroupInfo>& ffInfos) {
            // Group FFs by quantized depth
            std::map<int, std::vector<ClockInst*>> depthGroups;
            for (auto& fi : ffInfos) {
              depthGroups[fi.depth].push_back(fi.inst);
            }

            for (auto& [depth, ffs] : depthGroups) {
              if (depth == 0) {
                // Direct connection to driver subnet
                for (auto* ff : ffs) {
                  driverSubNet->addInst(*ff);
                  ++branchDirectFFs;
                }
              } else {
                // Determine group tier from majority vote of member FFs
                int grpBottom = 0, grpUpper = 0;
                for (auto* ff : ffs) {
                  int ft = cts3dDb_->getInstTier(ff->getDbInst());
                  if (ft == 0) ++grpBottom; else ++grpUpper;
                }
                const int grpTier = (grpBottom > grpUpper) ? 0 : 1;

                // Filter out cross-tier minority FFs before building the chain.
                // Minority FFs with a different tier than grpTier are direct-connected
                // to driverSubNet (the leaf net) instead of the group chain output net.
                // Rationale: a small chain-output net spanning two tiers triggers
                // GRT-0183 (heap underflow) because local HB via capacity is exhausted.
                // The leaf net (larger, many sinks spread across the die) handles
                // cross-tier routing flexibly, exactly as balanced CTS already does.
                // No extra chain buffers are added for minority FFs → buffer count
                // stays the same, keeping the LP problem efficient.
                std::vector<ClockInst*> sameTierFFs;
                for (auto* ff : ffs) {
                  int ft = cts3dDb_->getInstTier(ff->getDbInst());
                  if (ft < 0) ft = grpTier;  // unknown tier → assume same as group
                  if (ft != grpTier) {
                    // Cross-tier minority FF: return to leaf net (original placement)
                    driverSubNet->addInst(*ff);
                    ++branchDirectFFs;
                  } else {
                    sameTierFFs.push_back(ff);
                  }
                }

                if (sameTierFFs.empty()) {
                  // All FFs in this depth bucket were cross-tier; nothing to chain.
                  totalBranchSinks_debug += ffs.size();
                  continue;
                }

                // Compute centroid of same-tier FFs only (accurate chain placement)
                long long sumX = 0, sumY = 0;
                for (auto* ff : sameTierFFs) {
                  sumX += ff->getX();
                  sumY += ff->getY();
                }
                const int cx = static_cast<int>(sumX / sameTierFFs.size());
                const int cy = static_cast<int>(sumY / sameTierFFs.size());
                const std::string grpBufMaster = cts3dDb_->getBufferForTier(
                    options_->getRootBuffer(), grpTier);

                // Build shared chain: driver → grp_buf_0 → ... → grp_buf_{depth-1} → FFs
                // All chain buffers and sink FFs are the same tier → no cross-tier
                // chain nets → no GRT-0183.
                // V59: Sequential fixed-spacing placement along driver→centroid vector.
                // minStep = bufWidth → buffer bodies cannot overlap.
                // Replaces equal-fraction interpolation which caused overlap when
                // span < depth × bufWidth (DRT unroutable, GUI-0070).
                int grpBufW = 0;
                {
                  odb::dbMaster* m = db_->findMaster(grpBufMaster.c_str());
                  if (m) grpBufW = m->getWidth();
                }
                int grpSpanX = std::abs(cx - driverInst->getX());
                int grpSpanY = std::abs(cy - driverInst->getY());
                int grpSpan = std::max(grpSpanX, grpSpanY);
                int grpMinStep = std::max(grpBufW, 1);
                double grpFracStep = (grpSpan > 0)
                    ? static_cast<double>(grpMinStep) / grpSpan : 1.0;

                ClockSubNet* cur = driverSubNet;
                int prevGrpBufX = driverInst->getX();
                for (int d = 0; d < depth; ++d) {
                  double frac = grpFracStep * (d + 1);
                  if (frac >= 1.0) break;  // no room → stop chain
                  int bx = driverInst->getX()
                      + static_cast<int>((cx - driverInst->getX()) * frac);
                  int by = driverInst->getY()
                      + static_cast<int>((cy - driverInst->getY()) * frac);
                  {
                    Point<double> rawLoc(static_cast<double>(bx) / wireSegmentUnit_,
                                        static_cast<double>(by) / wireSegmentUnit_);
                    Point<double> legalLoc = legalizeOneBuffer(rawLoc, grpBufMaster);
                    bx = static_cast<int>(legalLoc.getX() * wireSegmentUnit_);
                    by = static_cast<int>(legalLoc.getY() * wireSegmentUnit_);
                  }
                  // Post-snap X-only overlap check: legalizeOneBuffer can snap to
                  // a different ROW, breaking the pre-snap spacing guarantee.
                  if (grpBufW > 0 && std::abs(bx - prevGrpBufX) <= grpBufW) {
                    break;  // overlap → stop chain
                  }
                  prevGrpBufX = bx;

                  // use myGrpId (global counter) instead of idx
                  // to avoid collision across sub-trees.
                  const std::string bName = tierPrefix_ + "clkbuf_grp_"
                      + std::to_string(myGrpId) + "_c"
                      + std::to_string(numSinks) + "_d"
                      + std::to_string(depth) + "_"
                      + std::to_string(d);
                  ClockInst& grpBuf = clock_.addClockBuffer(
                      bName, grpBufMaster, bx, by);
                  cts3dDb_->setClockInstTier(grpBuf, grpTier);
                  addTreeLevelBuffer(&grpBuf);

                  cur->addInst(grpBuf);
                  ClockSubNet& next = clock_.addSubNet(
                      tierPrefix_ + "clknet_grp_" + std::to_string(myGrpId) + "_c"
                      + std::to_string(numSinks) + "_d"
                      + std::to_string(depth) + "_"
                      + std::to_string(d));
                  next.addInst(grpBuf);
                  cur = &next;
                  ++branchGroupBufs;
                }

                // Connect same-tier FFs to final subnet (guaranteed single-tier net)
                cur->setLeafLevel(true);
                for (auto* ff : sameTierFFs) {
                  cur->addInst(*ff);
                  ++branchGroupedFFs;
                }
              }

              // Debug stats
              totalLpDelay_debug += 0;  // grouped: tracked at group level
              totalBranchSinks_debug += ffs.size();
            }
          };

          // GH-Tree Cascaded Delay Tap.
          // Replaces parallel grouped delay chains with cascaded tap structure.
          // Parallel chains (V49/V50):
          //   leaf -> [chain_d1]           -> d1_FFs    (1 buf, 1 net)
          //   leaf -> [buf -> buf]         -> d2_FFs    (2 bufs, 2 nets)
          //   leaf -> [buf -> buf -> buf]  -> d3_FFs    (3 bufs, 3 nets)
          //   Total: 6 bufs, 6 nets for maxDepth=3
          //
          // Cascaded tap (V52 GH-Tree):
          //   leaf -> tap1 -> tap2 -> tap3
          //             |       |       |
          //          d1_FFs  d2_FFs  d3_FFs
          //   Total: 3 bufs, 3 nets for maxDepth=3  (50% reduction)
          //
          // Same depth delivery: depth=k FF passes through exactly k tap buffers.
          // d0_FFs connect directly to driverSubNet (no tap buffer).
          // Same-tier purification (V49a) applied before cascaded chain.
          auto buildGHSubTree = [&](ClockInst* driverInst,
                                     ClockSubNet* driverSubNet,
                                     std::vector<FfGroupInfo>& ffInfos) {
            // Group FFs by quantized depth
            std::map<int, std::vector<ClockInst*>> depthGroups;
            for (auto& fi : ffInfos) {
              depthGroups[fi.depth].push_back(fi.inst);
            }

            // depth=0 FFs: direct connection to driver subnet
            if (depthGroups.count(0)) {
              for (auto* ff : depthGroups[0]) {
                driverSubNet->addInst(*ff);
                ++branchDirectFFs;
              }
            }

            // Find max depth with non-empty FFs
            int maxDepth = 0;
            for (auto& [d, ffs] : depthGroups) {
              if (d > 0 && !ffs.empty()) maxDepth = d;
            }
            if (maxDepth == 0) return;  // All depth=0, no sub-tree needed

            // Determine chain tier from ALL depth>0 FFs (majority vote)
            int chainBottom = 0, chainUpper = 0;
            for (auto& [depth, ffs] : depthGroups) {
              if (depth == 0) continue;
              for (auto* ff : ffs) {
                int ft = cts3dDb_->getInstTier(ff->getDbInst());
                if (ft == 0) ++chainBottom; else ++chainUpper;
              }
            }
            const int chainTier = (chainBottom > chainUpper) ? 0 : 1;

            // same-tier purification: cross-tier minority FFs -> depth=0
            for (auto& [depth, ffs] : depthGroups) {
              if (depth == 0) continue;
              auto it = ffs.begin();
              while (it != ffs.end()) {
                int ft = cts3dDb_->getInstTier((*it)->getDbInst());
                if (ft < 0) ft = chainTier;  // unknown -> assume same
                if (ft != chainTier) {
                  driverSubNet->addInst(**it);
                  ++branchDirectFFs;
                  it = ffs.erase(it);
                } else {
                  ++it;
                }
              }
            }

            // Recompute maxDepth after purification
            maxDepth = 0;
            for (auto& [d, ffs] : depthGroups) {
              if (d > 0 && !ffs.empty()) maxDepth = d;
            }
            if (maxDepth == 0) return;

            // Cost-aware adaptive depth selection
            // Scan from maxDepth downward. For each depth gap (empty taps),
            // check if the FFs beyond the gap justify the buffer cost.
            // If cost > benefit, reassign those FFs to the lower depth.
            // Proposition: reassign when ff_count * buf_delay * gap < gap * costRatio
            //   i.e., ff_count * buf_delay < costRatio  (per-gap-tap basis)
            {
              // Collect occupied depths in ascending order
              std::vector<int> occupiedDepths;
              for (auto& [d, ffs] : depthGroups) {
                if (d > 0 && !ffs.empty()) {
                  occupiedDepths.push_back(d);
                }
              }
              // Scan from highest depth downward
              bool pruned = true;
              while (pruned && occupiedDepths.size() > 1) {
                pruned = false;
                int highD = occupiedDepths.back();
                int lowD  = occupiedDepths[occupiedDepths.size() - 2];
                int gap   = highD - lowD;  // empty taps between lowD and highD
                int ffCount = static_cast<int>(depthGroups[highD].size());

                // Cost:    gap empty tap buffers (one buf+net per empty depth)
                // Benefit: ffCount FFs get (highD - lowD) * buf_delay more delivery
                // Reassign if ffCount is too small to justify the gap cost.
                // Threshold: ffCount * perFfBufDelay_ (ns) < adaptiveDepthCostRatio_ * gap_cost
                // Simplified: ffCount < costRatio * gap  (1 FF per empty tap)
                bool shouldReassign = false;
                if (gap > 1) {
                  // Gap has (gap-1) empty taps. Reassign if too few FFs.
                  int emptyTaps = gap - 1;
                  shouldReassign = (ffCount <= emptyTaps * adaptiveDepthCostRatio_);
                }
                // Also reassign if only 1 FF at the highest depth and gap >= 2
                if (ffCount == 1 && gap >= 2) {
                  shouldReassign = true;
                }

                if (shouldReassign) {
                  // Move highD FFs to lowD
                  for (auto* ff : depthGroups[highD]) {
                    depthGroups[lowD].push_back(ff);
                  }
                  debugPrint(logger_, CTS, "HTree", 1,
                      "Adaptive depth: reassign {} FFs depth {} -> {} "
                      "(gap={}, saved {} empty taps)",
                      ffCount, highD, lowD, gap, gap - 1);
                  depthGroups.erase(highD);
                  occupiedDepths.pop_back();
                  pruned = true;
                }
              }
              // Recompute maxDepth after adaptive pruning
              maxDepth = 0;
              for (auto& [d, ffs] : depthGroups) {
                if (d > 0 && !ffs.empty()) maxDepth = d;
              }
              if (maxDepth == 0) return;
            }

            // Compute centroid of ALL depth>0 same-tier FFs
            long long sumX = 0, sumY = 0;
            int nDownstream = 0;
            for (auto& [depth, ffs] : depthGroups) {
              if (depth == 0) continue;
              for (auto* ff : ffs) {
                sumX += ff->getX();
                sumY += ff->getY();
                ++nDownstream;
              }
            }
            if (nDownstream == 0) return;
            const int cx = static_cast<int>(sumX / nDownstream);
            const int cy = static_cast<int>(sumY / nDownstream);
            const std::string chainBufMaster = cts3dDb_->getBufferForTier(
                options_->getRootBuffer(), chainTier);

            // Build cascaded delay tap: driver -> tap1 -> tap2 -> ... -> tap_maxDepth
            // V59: Sequential fixed-spacing placement along driver→centroid vector.
            // minStep = bufWidth → buffer bodies cannot overlap.
            // Replaces physMax pre-cap + equal-fraction interpolation + post-snap check.
            int tapBufW = 0;
            {
              odb::dbMaster* m = db_->findMaster(chainBufMaster.c_str());
              if (m) tapBufW = m->getWidth();
            }
            int tapSpanX = std::abs(cx - driverInst->getX());
            int tapSpanY = std::abs(cy - driverInst->getY());
            int tapSpan = std::max(tapSpanX, tapSpanY);
            int tapMinStep = std::max(tapBufW, 1);
            double tapFracStep = (tapSpan > 0)
                ? static_cast<double>(tapMinStep) / tapSpan : 1.0;

            // Check: if even 1 buffer doesn't fit, skip chain entirely.
            if (tapFracStep >= 1.0) {
              logger_->info(CTS, 581,
                  "Cascaded TAP: chain skipped (span={}DBU, "
                  "bufWidth={}DBU — no room for tap buffers)",
                  tapSpan, tapBufW);
              return;
            }

            ClockSubNet* curNet = driverSubNet;
            int prevTapX = driverInst->getX();
            for (int d = 1; d <= maxDepth; ++d) {
              double frac = tapFracStep * d;
              if (frac >= 1.0) {
                // No room for more buffers — reassign remaining FFs to current depth
                for (int rd = d; rd <= maxDepth; ++rd) {
                  if (depthGroups.count(rd)) {
                    for (auto* ff : depthGroups[rd]) {
                      curNet->addInst(*ff);
                      ++branchGroupedFFs;
                    }
                  }
                }
                maxDepth = d - 1;
                break;
              }
              int bx = driverInst->getX()
                  + static_cast<int>((cx - driverInst->getX()) * frac);
              int by = driverInst->getY()
                  + static_cast<int>((cy - driverInst->getY()) * frac);

              // Grid-snap via legalizeOneBuffer()
              {
                Point<double> rawLoc(static_cast<double>(bx) / wireSegmentUnit_,
                                     static_cast<double>(by) / wireSegmentUnit_);
                Point<double> legalLoc = legalizeOneBuffer(rawLoc, chainBufMaster);
                bx = static_cast<int>(legalLoc.getX() * wireSegmentUnit_);
                by = static_cast<int>(legalLoc.getY() * wireSegmentUnit_);
              }

              // Post-snap X-only overlap check
              if (tapBufW > 0 && std::abs(bx - prevTapX) <= tapBufW) {
                // Overlap — reassign remaining FFs and stop
                for (int rd = d; rd <= maxDepth; ++rd) {
                  if (depthGroups.count(rd)) {
                    for (auto* ff : depthGroups[rd]) {
                      curNet->addInst(*ff);
                      ++branchGroupedFFs;
                    }
                  }
                }
                maxDepth = d - 1;
                break;
              }
              prevTapX = bx;

              // Create tap buffer
              const std::string bName = tierPrefix_ + "clkbuf_ghtap_"
                  + std::to_string(myGrpId) + "_c"
                  + std::to_string(numSinks) + "_t"
                  + std::to_string(d);
              ClockInst& tapBuf = clock_.addClockBuffer(
                  bName, chainBufMaster, bx, by);
              cts3dDb_->setClockInstTier(tapBuf, chainTier);
              addTreeLevelBuffer(&tapBuf);

              // Connect tap buffer input to previous net
              curNet->addInst(tapBuf);

              // Create tap output net
              ClockSubNet& tapNet = clock_.addSubNet(
                  tierPrefix_ + "clknet_ghtap_" + std::to_string(myGrpId) + "_c"
                  + std::to_string(numSinks) + "_t"
                  + std::to_string(d));
              tapNet.addInst(tapBuf);

              // Connect depth=d FFs to this tap net
              if (depthGroups.count(d) && !depthGroups[d].empty()) {
                for (auto* ff : depthGroups[d]) {
                  tapNet.addInst(*ff);
                  ++branchGroupedFFs;
                }
              }

              // Mark final tap as leaf level
              if (d == maxDepth) {
                tapNet.setLeafLevel(true);
              }

              curNet = &tapNet;
              ++branchGroupBufs;
            }

            totalBranchSinks_debug += nDownstream;
          };

          // Iterate sinks in this branch (same pattern as per-FF relay)
          const std::vector<Point<double>>& sinkLocs
              = leafTopology.getBranchSinksLocations(idx);
          for (const auto& loc : sinkLocs) {
            // Epsilon-tolerant sink lookup (V37 CTS-0080 fix)
            // FloatFix: eps=1e-4 (float→double mismatch reaches ~1e-5 for large coords)
            auto sinkIt = mapLocationToSink_.find(loc);
            if (sinkIt == mapLocationToSink_.end()) {
              constexpr double eps = 0.01;
              auto hint = mapLocationToSink_.lower_bound(
                  Point<double>(loc.getX() - eps, loc.getY() - eps));
              for (auto it = hint; it != mapLocationToSink_.end(); ++it) {
                if (it->first.getX() > loc.getX() + eps) break;
                if (std::abs(it->first.getX() - loc.getX()) < eps
                    && std::abs(it->first.getY() - loc.getY()) < eps) {
                  sinkIt = it;
                  break;
                }
              }
            }
            if (sinkIt == mapLocationToSink_.end()) {
              // Second-level clustering or legalizeOneBuffer can shift centroids
              // beyond eps=1e-4, causing lookup failure. Skip this sink —
              // it connects directly to the leaf net (depth=0, no grouped delay).
              debugPrint(logger_, CTS, "HTree", 1,
                  "Grouped: Sink not found at ({:.10f}, {:.10f}), skipping "
                  "(connects to leaf net directly)",
                  loc.getX(), loc.getY());
              continue;
            }
            ClockInst* sinkInst = sinkIt->second;

            // Check if cluster buffer or individual FF
            auto cIt = clusterBufTarget_.find(sinkInst);
            if (cIt != clusterBufTarget_.end()) {
              // Cluster buffer: add to branch subnet, then group its member FFs
              subNet->addInst(*sinkInst);

              ClockSubNet* leafSubNet = nullptr;
              clock_.forEachSubNet([&](ClockSubNet& sn) {
                if (sn.getDriver() == sinkInst) {
                  leafSubNet = &sn;
                }
              });

              if (leafSubNet != nullptr) {
                std::vector<ClockInst*> memberFFs;
                leafSubNet->forEachSink([&](ClockInst* ff) {
                  memberFFs.push_back(ff);
                });
                // Remove all sinks (will re-add via group chains)
                std::set<ClockInst*> toRemove(
                    memberFFs.begin(), memberFFs.end());
                leafSubNet->removeSinks(toRemove);

                // Collect FF info and build group chains
                std::vector<FfGroupInfo> ffInfos;
                collectFFs(sinkInst, leafSubNet, ffInfos, memberFFs);
                applyHoldBudgetCap(ffInfos);
                // GH-Tree cascaded tap or legacy parallel chains
                if (enableGhTree_) {
                  buildGHSubTree(sinkInst, leafSubNet, ffInfos);
                } else {
                  buildGroupChains(sinkInst, leafSubNet, ffInfos);
                }
              } else {
                ++branchDirectFFs;
              }
            } else {
              // Individual FF (no clustering): treat as single-FF group
              std::string ffName = getFFName(sinkInst);
              double ffTarget = cts3dDb_->getSkewTarget(ffName);
              // Debug: trace individual FF target lookup
              if (ffTarget > 0.01) {
                debugPrint(logger_, CTS, "HTree", 1,
                    "collectFFs individual: '{}' target={:.4f}ns",
                    ffName, ffTarget);
              }
              // globalMean fix: absolute LP target
              double ffDelta = ffTarget;
              int depth = 0;
              if (ffDelta > groupedDelayMinDelta_) {
                depth = static_cast<int>(
                    std::round(ffDelta / perFfBufDelay_));
                depth = std::max(0, std::min(depth, groupedDelayMaxDepth_));
              }
              // Removed redundant holdSkip (same as above).
              // V58 Fix C: Macro sinks (SRAM) get depth=0.
              odb::dbInst* dInst = sinkInst->getDbInst();
              if (isMacroBlockInst(dInst)) {
                depth = 0;
              }
              if (depth == 0) {
                subNet->addInst(*sinkInst);
                ++branchDirectFFs;
              } else {
                // Single FF needs its own chain
                int ft = cts3dDb_->getInstTier(sinkInst->getDbInst());
                const int ffTier = (ft >= 0) ? ft : branch_tier;
                const std::string ffBufMaster = cts3dDb_->getBufferForTier(
                    options_->getRootBuffer(), ffTier);
                // V59: Sequential fixed-spacing placement along driver→FF vector.
                // minStep = bufWidth → buffer bodies cannot overlap.
                int bufW = 0;
                {
                  odb::dbMaster* m = db_->findMaster(ffBufMaster.c_str());
                  if (m) bufW = m->getWidth();
                }
                int ffSpanX = std::abs(sinkInst->getX() - subNet->getDriver()->getX());
                int ffSpanY = std::abs(sinkInst->getY() - subNet->getDriver()->getY());
                int ffSpan = std::max(ffSpanX, ffSpanY);
                int ffMinStep = std::max(bufW, 1);
                double ffFracStep = (ffSpan > 0)
                    ? static_cast<double>(ffMinStep) / ffSpan : 1.0;

                if (ffFracStep >= 1.0 || depth == 0) {
                  // No room for any buffer → connect directly
                  subNet->addInst(*sinkInst);
                  ++branchDirectFFs;
                  ++numSinks;
                  continue;
                }
                ClockSubNet* cur = subNet;
                int prevFfBufX = subNet->getDriver()->getX();
                for (int d = 0; d < depth; ++d) {
                  double frac = ffFracStep * (d + 1);
                  if (frac >= 1.0) {
                    // No room → connect FF to current net and stop
                    cur->addInst(*sinkInst);
                    cur->setLeafLevel(true);
                    depth = d;
                    break;
                  }
                  int bx = subNet->getDriver()->getX()
                      + static_cast<int>(
                          (sinkInst->getX() - subNet->getDriver()->getX()) * frac);
                  int by = subNet->getDriver()->getY()
                      + static_cast<int>(
                          (sinkInst->getY() - subNet->getDriver()->getY()) * frac);
                  // Grid-snap per-FF grouped delay buffer
                  {
                    Point<double> rawLoc(static_cast<double>(bx) / wireSegmentUnit_,
                                        static_cast<double>(by) / wireSegmentUnit_);
                    Point<double> legalLoc = legalizeOneBuffer(rawLoc, ffBufMaster);
                    bx = static_cast<int>(legalLoc.getX() * wireSegmentUnit_);
                    by = static_cast<int>(legalLoc.getY() * wireSegmentUnit_);
                  }
                  // Post-snap X-only overlap check
                  if (bufW > 0 && std::abs(bx - prevFfBufX) <= bufW) {
                    cur->addInst(*sinkInst);
                    cur->setLeafLevel(true);
                    depth = d;
                    break;
                  }
                  prevFfBufX = bx;
                  const std::string bName = tierPrefix_ + "clkbuf_grp_"
                      + std::to_string(myGrpId) + "_s"
                      + std::to_string(numSinks) + "_"
                      + std::to_string(d);
                  ClockInst& grpBuf = clock_.addClockBuffer(
                      bName, ffBufMaster, bx, by);
                  cts3dDb_->setClockInstTier(grpBuf, ffTier);
                  addTreeLevelBuffer(&grpBuf);
                  cur->addInst(grpBuf);
                  ClockSubNet& next = clock_.addSubNet(
                      tierPrefix_ + "clknet_grp_" + std::to_string(myGrpId) + "_s"
                      + std::to_string(numSinks) + "_"
                      + std::to_string(d));
                  next.addInst(grpBuf);
                  cur = &next;
                  ++branchGroupBufs;
                }
                cur->setLeafLevel(true);
                cur->addInst(*sinkInst);
                ++branchGroupedFFs;
              }
            }
            ++numSinks;
          }

          numDelayBufs += branchGroupBufs;
          ++totalBranches_debug;
          if (branchGroupedFFs > 0) ++firedBranches_debug;
          totalDelivered_debug += branchGroupBufs * perFfBufDelay_;
          debugPrint(logger_, CTS, "HTree", 1,
              "{} branch {:2d}: grouped={} FFs ({} bufs), "
              "direct={}, holdSkip={}, branchMean={:.4f}ns",
              enableGhTree_ ? "ghtree" : "grouped", idx,
              branchGroupedFFs, branchGroupBufs,
              branchDirectFFs, branchHoldSkip, absTgt);
          return;  // skip legacy per-branch N-buffer code
        }

        // Per-FF relay buffer path (rewrote V39).
        // BUG #1 fix: V39 operated at cluster-buffer granularity (mean target),
        //   now iterates member FFs individually for per-FF relay decisions.
        // BUG #2 fix: V39 placed relays at FF position (unplanned wire delay),
        //   now places relays at interpolated positions between leaf buffer and FF.
        // BUG #3 fix: Hold guard now checks member FF hold budgets (min).
        //
        // Structure: branch → leaf_buffer → [relay_0] → [relay_1] → FF
        //            branch → leaf_buffer → FF  (zero target, direct)
        if (enablePerFfRelay_ && cts3dDb_ != nullptr
            && cts3dDb_->hasSkewTargets() && perFfBufDelay_ > 1e-9) {
          subNet->setLeafLevel(true);
          int branchRelayFFs = 0, branchRelayBufs = 0;
          int branchDirectFFs = 0, branchHoldSkip = 0;

          const int branch_tier = cts3dDb_->getDominantTier(
              leafTopology.getBranchSinksLocations(idx), mapLocationToSink_);
          const std::string bufMaster = cts3dDb_->getBufferForTier(
              options_->getRootBuffer(), branch_tier);

          // Helper lambda to process a single FF with relay decision
          auto processOneFf = [&](ClockInst* ffInst, ClockInst* driverInst,
                                  ClockSubNet* driverSubNet) {
            // Get per-FF target (NOT cluster mean)
            const double ffTarget = [&]() -> double {
              if (cts3dDb_ == nullptr || !cts3dDb_->hasSkewTargets()) return 0.0;
              std::string nm = ffInst->getName();
              const auto sp = nm.rfind('/');
              if (sp != std::string::npos) nm = nm.substr(0, sp);
              return cts3dDb_->getSkewTarget(nm);
            }();
            // globalMean fix: absolute LP target
            const double ffDelta = ffTarget;
            // Pre-compute dx/dy for Elmore model (moved before nRelay decision)
            const int dx = ffInst->getX() - driverInst->getX();
            const int dy = ffInst->getY() - driverInst->getY();

            // Elmore x_useful relay positioning.
            // When wire RC params are set (CTS_RELAY_RW_PER_UM / CTS_RELAY_CW_PER_UM),
            // compute exact relay position so 1 buffer delivers target_delta:
            //   delta(x) = d_buf - rc * x * (L - x)   [ps]
            //   x_useful = [L - sqrt(L^2 - 4*(d_buf-target)/rc)] / 2  [um]
            // where rc = rw[kOhm/um] * cw[fF/um] = ps/um^2 (kOhm*fF = 1ps).
            // For target > d_buf: fall back to multi-relay equal spacing.
            // For rc == 0: fall back to old integer-quantized V41 method.
            int nRelay = 0;
            double relayFrac = 0.5;  // position fraction for 1-relay Elmore case
            if (ffDelta > perFfMinTarget_) {
              const double d_buf_ps  = perFfBufDelay_ * 1000.0;  // ns -> ps
              const double target_ps = ffDelta * 1000.0;          // ns -> ps
              const double rc = relayRwKOhmPerUm_ * relayCwFfPerUm_;  // ps/um^2
              const bool elmore = (rc > 1e-12
                  && relayDbuPerUm_ > 1e-6
                  && (std::abs(dx) + std::abs(dy)) > 0);
              if (elmore) {
                const double L_um = (std::abs(dx) + std::abs(dy))
                    / relayDbuPerUm_;
                if (target_ps <= d_buf_ps && L_um > 1e-6) {
                  // 1 relay at x_useful covers target exactly
                  const double disc = L_um * L_um
                      - 4.0 * (d_buf_ps - target_ps) / rc;
                  if (disc >= 0.0) {
                    const double x_useful = (L_um - std::sqrt(disc)) / 2.0;
                    relayFrac = std::max(0.05,
                        std::min(0.95, x_useful / L_um));
                  }
                  // disc < 0 means target > max wire contribution → midpoint
                  nRelay = 1;
                } else if (L_um > 1e-6) {
                  // target > d_buf: multiple relays at equal spacing
                  nRelay = std::max(1, std::min(perFfMaxRelay_,
                      static_cast<int>(
                          std::round(target_ps / d_buf_ps))));
                }
              } else {
                // Fallback: integer quantization (V41 method)
                nRelay = static_cast<int>(
                    std::round(ffDelta / perFfBufDelay_));
                nRelay = std::max(0, std::min(nRelay, perFfMaxRelay_));
              }
            }

            // Removed redundant holdSkip (same as above).

            // Debug stats (use delta, not absolute, for accurate coverage)
            totalLpDelay_debug += std::max(0.0, ffDelta);
            ++totalBranchSinks_debug;

            if (nRelay > 0) {
              // dx/dy pre-computed above (before nRelay decision).
              // Single relay uses Elmore relayFrac; multi-relay uses equal spacing.
              // Get FF's tier for relay tier assignment.
              // Relay must be on the SAME tier as its target FF to avoid cross-tier
              // relay→FF nets. Cross-tier relay→FF nets cause GRT-0183 heap underflow
              // in 3D maze routing (observed in V43/V46/V47 for minority-tier FFs in
              // upper-dominant clusters). Using branch_tier (cluster dominant tier)
              // was wrong when FF is on the minority tier of the cluster.
              const int ff_tier_relay = cts3dDb_->getInstTier(ffInst->getDbInst());
              const int relay_tier = (ff_tier_relay >= 0) ? ff_tier_relay : branch_tier;

              ClockSubNet* cur = driverSubNet;
              for (int r = 0; r < nRelay; ++r) {
                // 1-relay → Elmore x_useful fraction; else equal spacing
                const double frac = (nRelay == 1)
                    ? relayFrac
                    : static_cast<double>(r + 1) / (nRelay + 1);
                const int rx = driverInst->getX()
                    + static_cast<int>(dx * frac);
                const int ry = driverInst->getY()
                    + static_cast<int>(dy * frac);

                const std::string rName = tierPrefix_ + "clkbuf_relay_"
                    + std::to_string(idx) + "_"
                    + std::to_string(numSinks) + "_"
                    + std::to_string(r);
                ClockInst& relay = clock_.addClockBuffer(
                    rName, bufMaster, rx, ry);
                cts3dDb_->setClockInstTier(relay, relay_tier);
                addTreeLevelBuffer(&relay);
                cur->addInst(relay);
                ClockSubNet& next = clock_.addSubNet(
                    tierPrefix_ + "clknet_relay_" + std::to_string(idx) + "_"
                    + std::to_string(numSinks) + "_"
                    + std::to_string(r));
                next.addInst(relay);
                cur = &next;
              }
              cur->setLeafLevel(true);
              cur->addInst(*ffInst);
              ++branchRelayFFs;
              branchRelayBufs += nRelay;
              totalDelivered_debug += nRelay * perFfBufDelay_;
            } else {
              driverSubNet->addInst(*ffInst);
              ++branchDirectFFs;
            }
            ++numSinks;
          };

          const std::vector<Point<double>>& sinkLocs
              = leafTopology.getBranchSinksLocations(idx);
          for (const auto& loc : sinkLocs) {
            // Epsilon-tolerant sink lookup (V37 CTS-0080 fix reused)
            // FloatFix: eps=1e-4 (float→double mismatch reaches ~1e-5 for large coords)
            auto sinkIt = mapLocationToSink_.find(loc);
            if (sinkIt == mapLocationToSink_.end()) {
              constexpr double eps = 0.01;
              auto hint = mapLocationToSink_.lower_bound(
                  Point<double>(loc.getX() - eps, loc.getY() - eps));
              for (auto it = hint; it != mapLocationToSink_.end(); ++it) {
                if (it->first.getX() > loc.getX() + eps) break;
                if (std::abs(it->first.getX() - loc.getX()) < eps
                    && std::abs(it->first.getY() - loc.getY()) < eps) {
                  sinkIt = it;
                  break;
                }
              }
            }
            if (sinkIt == mapLocationToSink_.end()) {
              logger_->error(CTS, 425,
                  "Relay: Sink not found at ({:.10f}, {:.10f})",
                  loc.getX(), loc.getY());
            }
            ClockInst* sinkInst = sinkIt->second;

            // Check if this is a cluster buffer
            auto cIt = clusterBufTarget_.find(sinkInst);
            if (cIt != clusterBufTarget_.end()) {
              // Cluster buffer: add to branch subnet, then process
              // each member FF individually with per-FF relay chains.
              subNet->addInst(*sinkInst);

              // Find the cluster's leaf subnet
              ClockSubNet* leafSubNet = nullptr;
              clock_.forEachSubNet([&](ClockSubNet& sn) {
                if (sn.getDriver() == sinkInst) {
                  leafSubNet = &sn;
                }
              });

              if (leafSubNet != nullptr) {
                // Collect member FFs (can't modify while iterating)
                std::vector<ClockInst*> memberFFs;
                leafSubNet->forEachSink([&](ClockInst* ff) {
                  memberFFs.push_back(ff);
                });
                // Remove all sinks from leaf subnet (will re-add or relay)
                std::set<ClockInst*> toRemove(
                    memberFFs.begin(), memberFFs.end());
                leafSubNet->removeSinks(toRemove);

                // Process each member FF individually
                for (ClockInst* memberFF : memberFFs) {
                  processOneFf(memberFF, sinkInst, leafSubNet);
                }
              } else {
                // Leaf subnet not found (shouldn't happen)
                ++branchDirectFFs;
                ++numSinks;
              }
            } else {
              // Individual FF (no clustering): process directly
              processOneFf(sinkInst, sinkInst, subNet);
            }
          }

          numDelayBufs += branchRelayBufs;
          ++totalBranches_debug;
          if (branchRelayFFs > 0) ++firedBranches_debug;
          logger_->info(CTS, 424,
              "Relay branch {:2d}: relay={} FFs ({} bufs), "
              "direct={}, holdSkip={}, branchMean={:.4f}ns",
              idx, branchRelayFFs, branchRelayBufs,
              branchDirectFFs, branchHoldSkip, absTgt);
          return;  // skip legacy per-branch N-buffer code
        }

        // globalMean fix: absolute LP target
        const double delta = absTgt;
        int num_bufs = 0;
        // V58 Fix F: Skip per-branch delay buffers for MacroTree.
        // Macro tree latency is already too high (1,095ps vs register ~462ps).
        // Adding delay buffers widens the gap → worsens hold on reg→macro edges.
        if (type_ != TreeType::MacroTree
            && delta > 0.0 && singleBufDelay_ > 1e-9) {
          num_bufs = static_cast<int>(std::round(delta / singleBufDelay_));
          num_bufs = std::max(0, std::min(num_bufs, maxLeafDelayBufs_));
        }

        // Debug: accumulate LP-CTS gap stats (CTS-408 summary after loop)
        totalLpDelay_debug   += absTgt * branchCnt;
        totalBranchSinks_debug += branchCnt;
        ++totalBranches_debug;

        if (num_bufs > 0) {
            totalDelivered_debug += num_bufs * singleBufDelay_ * branchCnt;
            ++firedBranches_debug;

            const int branch_tier = cts3dDb_->getDominantTier(
                leafTopology.getBranchSinksLocations(idx), mapLocationToSink_);
            const std::string delayBufMaster = cts3dDb_->getBufferForTier(
                options_->getRootBuffer(), branch_tier);

            const int bx = branchPoint.getX() * wireSegmentUnit_;
            const int by = branchPoint.getY() * wireSegmentUnit_;

            // Build N-buffer chain in series
            ClockSubNet* curSub = subNet;
            for (int b = 0; b < num_bufs; ++b) {
              const std::string bufName = tierPrefix_ + "clkbuf_delay_"
                                          + std::to_string(idx) + "_"
                                          + std::to_string(b);
              ClockInst& buf = clock_.addClockBuffer(bufName, delayBufMaster, bx, by);
              cts3dDb_->setClockInstTier(buf, branch_tier);
              addTreeLevelBuffer(&buf);

              // buf is SINK of curSub (curSub already has a driver)
              curSub->addInst(buf);

              // Create new subnet; first addInst sets buf as its driver
              ClockSubNet& nextSub = clock_.addSubNet(
                  tierPrefix_ + "clknet_delay_" + std::to_string(idx) + "_" + std::to_string(b));
              nextSub.addInst(buf);  // buf becomes driver of nextSub

              curSub = &nextSub;
            }

            // Last subnet is the leaf subnet; connect FFs as sinks
            curSub->setLeafLevel(true);
            const std::vector<Point<double>>& sinkLocs
                = leafTopology.getBranchSinksLocations(idx);
            for (const Point<double>& loc : sinkLocs) {
              auto sinkIt = mapLocationToSink_.find(loc);
              // epsilon fallback for float precision mismatch
              // FloatFix: eps=1e-4 (float→double mismatch reaches ~1e-5 for large coords)
              if (sinkIt == mapLocationToSink_.end()) {
                constexpr double eps = 0.01;
                auto hint = mapLocationToSink_.lower_bound(
                    Point<double>(loc.getX() - eps, loc.getY() - eps));
                for (auto it = hint; it != mapLocationToSink_.end(); ++it) {
                  if (it->first.getX() > loc.getX() + eps) break;
                  if (std::abs(it->first.getX() - loc.getX()) < eps
                      && std::abs(it->first.getY() - loc.getY()) < eps) {
                    sinkIt = it;
                    break;
                  }
                }
              }
              if (sinkIt == mapLocationToSink_.end()) {
                // V58: eps failed — find absolute nearest (same pattern as CTS-0079).
                // Cannot skip: sink would be unconnected → SIGSEGV in writeDummyLoadsToDb.
                double bestDist = std::numeric_limits<double>::max();
                for (auto& [key, sink] : mapLocationToSink_) {
                  double d = std::abs(key.getX() - loc.getX())
                           + std::abs(key.getY() - loc.getY());
                  if (d < bestDist) {
                    bestDist = d;
                    sinkIt = mapLocationToSink_.find(key);
                  }
                }
                logger_->warn(CTS, 404,
                    "Delay buf branch: sink at ({:.6f}, {:.6f}) — nearest match dist={:.6f}",
                    loc.getX(), loc.getY(), bestDist);
              }
              curSub->addInst(*sinkIt->second);
              ++numSinks;
            }

            ++numDelayBufs;
            logger_->info(CTS, 402,
                          "Delay buffers inserted at leaf branch {}: "
                          "absTarget={:.4f}ns, globalMean={:.4f}ns, delta={:.4f}ns, N={}, buffer={}",
                          idx, absTgt, globalMeanTarget_, delta, num_bufs, delayBufMaster);
            return;  // skip normal sink connection below
          }

        // Branch has LP target but delta ≤ 0 → skipped by delta trigger.
        // These branches are BELOW globalMean; N=0 even though LP assigned non-zero target.
        // This reveals LP targets that are "wasted" (below average, can't be implemented
        // with delta-based approach). Key gap source if most targets fall here.
        if (absTgt > 0.0) {
          logger_->info(CTS, 406,
                        "  Branch {:2d} skipped: absTgt={:.4f}ns <= globalMean={:.4f}ns "
                        "(delta={:.4f}ns -> N=0, {} sinks; "
                        "LP wanted {:.4f}ns/sink but CTS delivers 0ns)",
                        idx, absTgt, globalMeanTarget_, delta, branchCnt, absTgt);
        }

        subNet->setLeafLevel(true);

        const std::vector<Point<double>>& sinkLocs
            = leafTopology.getBranchSinksLocations(idx);
        for (const Point<double>& loc : sinkLocs) {
          auto sinkIt = mapLocationToSink_.find(loc);
          // float/double precision mismatch from multi-level
          // computeBranchSinks() double→float→double roundtrip.
          // If exact lookup fails, search for nearest entry within epsilon.
          // FloatFix: eps=1e-4 (float→double mismatch reaches ~1e-5 for large coords)
          if (sinkIt == mapLocationToSink_.end()) {
            constexpr double eps = 0.01;
            auto hint = mapLocationToSink_.lower_bound(
                Point<double>(loc.getX() - eps, loc.getY() - eps));
            for (auto it = hint; it != mapLocationToSink_.end(); ++it) {
              if (it->first.getX() > loc.getX() + eps) break;
              if (std::abs(it->first.getX() - loc.getX()) < eps
                  && std::abs(it->first.getY() - loc.getY()) < eps) {
                sinkIt = it;
                debugPrint(logger_, CTS, "clustering", 1,
                           "Leaf branch {}: epsilon match ({:.10f},{:.10f}) -> "
                           "({:.10f},{:.10f}) = {}",
                           idx, loc.getX(), loc.getY(),
                           it->first.getX(), it->first.getY(),
                           it->second->getName());
                break;
              }
            }
          }
          if (sinkIt == mapLocationToSink_.end()) {
            // V58: warn+skip instead of fatal error.
            // Multi-level clustering float precision drift can exceed eps=1e-4
            // for large designs (10K+ sinks, 7+ H-tree levels).
            // Skipped sink stays on parent leaf net (depth=0, no grouped delay).
            logger_->warn(CTS, 80,
                          "Sink not found at ({:.10f}, {:.10f}), mapSize={} — "
                          "skipping (connects to leaf net directly)",
                          loc.getX(), loc.getY(), mapLocationToSink_.size());
            continue;
          }

          subNet->addInst(*sinkIt->second);
          ++numSinks;
        }
      });

  logger_->info(CTS, 35, " Number of sinks covered: {}.", numSinks);
  if (numDelayBufs > 0) {
    logger_->info(CTS, 403,
                  " {} delay buffers inserted for per-branch skew targeting.",
                  numDelayBufs);
  }

  // LP-CTS translation efficiency summary.
  // coverage = totalDelivered / totalLpDelay * 100%
  //   ~0%   → N-buffer chain never fires (all branches below globalMean, or absTgt=0)
  //   ~50%  → half of LP target translated (rounding + below-mean branches)
  //   ~100% → N-buffer chain faithfully implements LP targets
  // This is the primary metric for diagnosing LP→CTS pipeline effectiveness.
  if (totalBranchSinks_debug > 0 && totalLpDelay_debug > 1e-12) {
    const double coverage = 100.0 * totalDelivered_debug / totalLpDelay_debug;
    logger_->info(CTS, 408,
                  "LP-CTS gap: LP_total={:.4f}ns ({:.4f}ns/sink), "
                  "delivered={:.4f}ns ({:.4f}ns/sink), coverage={:.1f}%, "
                  "fired={}/{} leaf branches",
                  totalLpDelay_debug,
                  totalLpDelay_debug / totalBranchSinks_debug,
                  totalDelivered_debug,
                  totalDelivered_debug / totalBranchSinks_debug,
                  coverage,
                  firedBranches_debug, totalBranches_debug);
  }

  // Per-FF relay summary log
  if (enablePerFfRelay_ && totalBranchSinks_debug > 0) {
    const double cov = totalLpDelay_debug > 1e-12
        ? 100.0 * totalDelivered_debug / totalLpDelay_debug : 0.0;
    logger_->info(CTS, 426,
        "Per-FF relay summary: {} relay bufs across {}/{} branches, "
        "{} sinks, coverage={:.1f}%",
        numDelayBufs, firedBranches_debug, totalBranches_debug,
        numSinks, cov);
  }
}

void HTreeBuilder::createSingleBufferClockNet()
{
  logger_->report(" Building single-buffer clock net.");

  Point<double> center = sinkRegion_.getCenter();
  Point<double> legalCenter
      = legalizeOneBuffer(center, options_->getRootBuffer());
  sinkRegion_.setCenter(legalCenter);
  commitMoveLoc(center, legalCenter);
  const int centerX = legalCenter.getX() * wireSegmentUnit_;
  const int centerY = legalCenter.getY() * wireSegmentUnit_;
  ClockInst& rootBuffer = clock_.addClockBuffer(
      tierPrefix_ + "clkbuf_0", options_->getRootBuffer(), centerX, centerY);

  // clang-format off
  if (center != legalCenter) {
    debugPrint(logger_, CTS, "legalizer", 2, "createSingleBufferClockNet "
	       "legalizeOneBuffer {}clkbuf_0: {} => {}", tierPrefix_, center, legalCenter);
  }
  // clang-format on

  addTreeLevelBuffer(&rootBuffer);
  ClockSubNet& clockSubNet = clock_.addSubNet(tierPrefix_ + "clknet_0");
  clockSubNet.setTreeLevel(0);
  clockSubNet.addInst(rootBuffer);

  clock_.forEachSink([&](ClockInst& inst) { clockSubNet.addInst(inst); });
}

void HTreeBuilder::plotSolution()
{
  auto name = std::string("plot_") + clock_.getName() + ".py";
  std::ofstream file(name);
  file << "import numpy as np\n";
  file << "import matplotlib.pyplot as plt\n";
  file << "import matplotlib.path as mpath\n";
  file << "import matplotlib.lines as mlines\n";
  file << "import matplotlib.patches as mpatches\n";
  file << "from matplotlib.collections import PatchCollection\n\n";

  clock_.forEachSink([&](const ClockInst& sink) {
    file << "plt.scatter(" << (double) sink.getX() / wireSegmentUnit_ << ", "
         << (double) sink.getY() / wireSegmentUnit_ << ", s=1)\n";
  });

  LevelTopology& topLevelTopology = topologyForEachLevel_.front();
  Point<double> topLevelBufferLoc = sinkRegion_.getCenter();
  topLevelTopology.forEachBranchingPoint(
      [&](unsigned idx, Point<double> branchPoint) {
        if (topLevelBufferLoc.getX() < branchPoint.getX()) {
          file << "plt.plot([" << topLevelBufferLoc.getX() << ", "
               << branchPoint.getX() << "], [" << topLevelBufferLoc.getY()
               << ", " << branchPoint.getY() << "], c = 'r')\n";
        } else {
          file << "plt.plot([" << branchPoint.getX() << ", "
               << topLevelBufferLoc.getX() << "], [" << branchPoint.getY()
               << ", " << topLevelBufferLoc.getY() << "], c = 'r')\n";
        }
      });

  for (int levelIdx = 1; levelIdx < topologyForEachLevel_.size(); ++levelIdx) {
    const LevelTopology& topology = topologyForEachLevel_[levelIdx];
    topology.forEachBranchingPoint([&](unsigned idx,
                                       Point<double> branchPoint) {
      unsigned parentIdx = topology.getBranchingPointParentIdx(idx);
      Point<double> parentPoint
          = topologyForEachLevel_[levelIdx - 1].getBranchingPoint(parentIdx);
      std::string color = "orange";
      if (levelIdx % 2 == 0) {
        color = "red";
      }

      if (parentPoint.getX() < branchPoint.getX()) {
        file << "plt.plot([" << parentPoint.getX() << ", " << branchPoint.getX()
             << "], [" << parentPoint.getY() << ", " << branchPoint.getY()
             << "], c = '" << color << "')\n";
      } else {
        file << "plt.plot([" << branchPoint.getX() << ", " << parentPoint.getX()
             << "], [" << branchPoint.getY() << ", " << parentPoint.getY()
             << "], c = '" << color << "')\n";
      }
    });
  }

  file << "plt.show()\n";
  file.close();
}

// print structures of Htree from top level buffer
// incluiding branch point locations, topology length and weighted sink lengths
void HTreeBuilder::printHTree()
{
  Point<double> topLevelBufferLoc = sinkRegion_.getCenter();
  logger_->report("HTree: top buf loc:{}", topLevelBufferLoc);
  for (int levelIdx = 0; levelIdx < topologyForEachLevel_.size(); ++levelIdx) {
    LevelTopology& topology = topologyForEachLevel_[levelIdx];

    for (unsigned idx = 0; idx < topology.getBranchingPointSize(); ++idx) {
      Point<double>& branchPoint = topology.getBranchingPoint(idx);
      unsigned parentIdx = topology.getBranchingPointParentIdx(idx);

      // clang-format off
      Point<double> parentPoint
          = (levelIdx == 0)
                ? topLevelBufferLoc
                : topologyForEachLevel_[levelIdx - 1].getBranchingPoint(
                      parentIdx);
      // clang-format on

      const std::vector<Point<double>>& sinks
          = topology.getBranchSinksLocations(idx);

      double leng = topology.getLength();
      // clang-format off
      logger_->report("HTree: level*{}* bufId*{}*: branchPt:{} topo len:{:0.3f}"
		      " dist to parent:{:0.3f} weighted sink len:{:0.3f} "
		      "parentPt:{}", levelIdx, idx, branchPoint, leng,
		      computeDist(branchPoint, parentPoint),
		      weightedDistance(branchPoint, branchPoint, sinks),
		      parentPoint);
      // clang-format on
    }
    logger_->report("-------------------------------------------------");
  }
}

SegmentBuilder::SegmentBuilder(const std::string& instPrefix,
                               const std::string& netPrefix,
                               const Point<double>& root,
                               const Point<double>& target,
                               const std::vector<unsigned>& techCharWires,
                               Clock& clock,
                               ClockSubNet& drivingSubNet,
                               const TechChar& techChar,
                               const unsigned techCharDistUnit,
                               TreeBuilder* tree,
                               odb::dbDatabase* db,
                               int targetTier)
    : instPrefix_(instPrefix),
      netPrefix_(netPrefix),
      root_(root),
      target_(target),
      techCharWires_(techCharWires),
      techChar_(&techChar),
      techCharDistUnit_(techCharDistUnit),
      clock_(&clock),
      drivingSubNet_(&drivingSubNet),
      tree_(tree),
      db_(db),
      targetTier_(targetTier)
{
}

void SegmentBuilder::build(const std::string& forceBuffer)
{
  const double lengthX = std::abs(root_.getX() - target_.getX());
  const bool isLowToHiX = root_.getX() < target_.getX();
  const bool isLowToHiY = root_.getY() < target_.getY();

  double connectionLength = 0.0;
  for (unsigned techCharWireIdx : techCharWires_) {
    const WireSegment& wireSegment = techChar_->getWireSegment(techCharWireIdx);
    const unsigned wireSegLen = wireSegment.getLength();
    for (int buffer = 0; buffer < wireSegment.getNumBuffers(); ++buffer) {
      const double location
          = wireSegment.getBufferLocation(buffer) * wireSegLen;
      connectionLength += location;

      double x = std::numeric_limits<double>::max();
      double y = std::numeric_limits<double>::max();
      if (connectionLength < lengthX) {
        y = root_.getY();
        x = (isLowToHiX) ? (root_.getX() + connectionLength)
                         : (root_.getX() - connectionLength);
      } else {
        x = target_.getX();
        y = (isLowToHiY) ? (root_.getY() + (connectionLength - lengthX))
                         : (root_.getY() - (connectionLength - lengthX));
      }

      const std::string buffMaster = !forceBuffer.empty()
                                         ? forceBuffer
                                         : wireSegment.getBufferMaster(buffer);
      // Replaced mapBufferMasterToTier with Cts3DDatabase
      const std::string tieredMaster
          = tree_->getCts3DDatabase()->getBufferForTier(buffMaster, targetTier_);
      Point<double> bufferLoc(x, y);
      Point<double> legalBufferLoc
          = tree_->legalizeOneBuffer(bufferLoc, tieredMaster);
      tree_->commitMoveLoc(bufferLoc, legalBufferLoc);
      ClockInst& newBuffer = clock_->addClockBuffer(
          instPrefix_ + std::to_string(numBufferLevels_),
          tieredMaster,
          legalBufferLoc.getX() * techCharDistUnit_,
          legalBufferLoc.getY() * techCharDistUnit_);
      // Tier assignment through SSOT (Cts3DDatabase)
      tree_->getCts3DDatabase()->setClockInstTier(newBuffer, targetTier_);

      // clang-format off
      if (bufferLoc != legalBufferLoc) {
	// adjust for cell movement
	connectionLength -= tree_->computeDist(bufferLoc, legalBufferLoc);
	debugPrint(getTree()->getLogger(), CTS, "legalizer", 2,
		   " SegmentBuilder::build {} TCId:{} bufId:{} connLen:{:0.1f}: "
		   "{} => {}", instPrefix_  + std::to_string(numBufferLevels_),
		   techCharWireIdx, buffer, connectionLength,
		   bufferLoc, legalBufferLoc);
      } else {
	debugPrint(getTree()->getLogger(), CTS, "legalizer", 2,
		   " SegmentBuilder::build {} TCId:{} bufId:{} connLen:{:0.1f}: "
		   "{}", instPrefix_  + std::to_string(numBufferLevels_),
		   techCharWireIdx, buffer, connectionLength, bufferLoc);
      }
      // clang-format on

      tree_->addTreeLevelBuffer(&newBuffer);

      drivingSubNet_->addInst(newBuffer);
      drivingSubNet_
          = &clock_->addSubNet(netPrefix_ + std::to_string(numBufferLevels_));
      drivingSubNet_->addInst(newBuffer);

      ++numBufferLevels_;
    }
    connectionLength += wireSegLen;
  }
}

void SegmentBuilder::forceBufferInSegment(const std::string& master)
{
  if (numBufferLevels_ != 0) {
    return;
  }

  // Replaced mapBufferMasterToTier with Cts3DDatabase
  const std::string tieredMaster
      = tree_->getCts3DDatabase()->getBufferForTier(master, targetTier_);
  ClockInst& newBuffer
      = clock_->addClockBuffer(instPrefix_ + "_f",
                               tieredMaster,
                               target_.getX() * techCharDistUnit_,
                               target_.getY() * techCharDistUnit_);
  // Tier assignment through SSOT (Cts3DDatabase)
  tree_->getCts3DDatabase()->setClockInstTier(newBuffer, targetTier_);
  tree_->addTreeLevelBuffer(&newBuffer);
  // clang-format off
  debugPrint(getTree()->getLogger(), CTS, "legalizer", 2,
	     "  forceBufferInSegment {}: {}", instPrefix_ + "_f", target_);
  // clang-format on

  drivingSubNet_->addInst(newBuffer);
  drivingSubNet_ = &clock_->addSubNet(netPrefix_ + "_leaf");
  drivingSubNet_->addInst(newBuffer);
  numBufferLevels_++;
}

// TARP: Timing-affinity clustering integration.
void HTreeBuilder::runTarpClustering(
    const std::vector<std::pair<float, float>>& points,
    const std::vector<const ClockInst*>& sinkInsts,
    unsigned hilbertBestGroupSize,
    float hilbertBestDiameter)
{
  const unsigned N = points.size();

  TarpParams tp;
  tp.alpha = tarpAlpha_;
  tp.d_buf = perFfBufDelay_;
  tp.max_depth = groupedDelayMaxDepth_;

  // Hilbert-calibrated cluster size: use grid search result as default.
  // Env var CTS_TARP_MAX_CLUSTER overrides if set.
  tp.max_cluster_size = (hilbertBestGroupSize > 0) ? hilbertBestGroupSize : 10;
  if (const char* e = std::getenv("CTS_TARP_MAX_CLUSTER"))
    tp.max_cluster_size = std::atoi(e);

  if (const char* e = std::getenv("CTS_TARP_TOP_K"))
    tp.top_k_edges = std::atoi(e);
  if (const char* e = std::getenv("CTS_TARP_SLACK_THRESH"))
    tp.slack_threshold_ratio = std::stod(e);
  if (const char* e = std::getenv("CTS_TARP_FM_PASSES"))
    tp.fm_max_passes = std::atoi(e);

  // Hilbert-calibrated diameter: use grid search result as default.
  // Env var CTS_TARP_MAX_DIAMETER_UM overrides if set.
  if (hilbertBestDiameter > 0) {
    tp.max_diameter = hilbertBestDiameter;  // already in wireSegmentUnit
  } else {
    tp.max_diameter = 50.0 * options_->getDbUnits() / wireSegmentUnit_;
  }
  if (const char* e = std::getenv("CTS_TARP_MAX_DIAMETER_UM"))
    tp.max_diameter = std::stod(e) * options_->getDbUnits() / wireSegmentUnit_;

  TarpClustering tarp(logger_, tp);

  // Add sinks with LP targets and hold budgets
  int hb_direct = 0, hb_downstream = 0, hb_fail_closed = 0;
  for (unsigned i = 0; i < N; ++i) {
    std::string nm = sinkInsts[i]->getName();
    const auto sp = nm.rfind('/');
    if (sp != std::string::npos) nm = nm.substr(0, sp);

    double lp_target = cts3dDb_->hasSkewTargets()
                           ? cts3dDb_->getSkewTarget(nm) : 0.0;
    double hold_budget = 0.0;
    auto it = perFfHoldBudget_.find(nm);
    if (it != perFfHoldBudget_.end()) {
      hold_budget = it->second;
      ++hb_direct;
    } else {
      // Clone/gate sink: trace downstream to find real FFs' hold budgets
      odb::dbInst* dbInst = sinkInsts[i]->getDbInst();
      bool found_downstream = false;
      if (dbInst) {
        for (auto* iterm : dbInst->getITerms()) {
          if (iterm->getIoType() != odb::dbIoType::OUTPUT) continue;
          odb::dbNet* outNet = iterm->getNet();
          if (!outNet) continue;
          for (auto* sink_iterm : outNet->getITerms()) {
            if (sink_iterm == iterm) continue;
            std::string ff_name = sink_iterm->getInst()->getName();
            auto ff_it = perFfHoldBudget_.find(ff_name);
            if (ff_it != perFfHoldBudget_.end()) {
              if (!found_downstream) {
                hold_budget = ff_it->second;
                found_downstream = true;
              } else {
                hold_budget = std::min(hold_budget, ff_it->second);
              }
            }
          }
        }
      }
      if (!found_downstream) {
        ++hb_fail_closed;
      } else {
        ++hb_downstream;
      }
    }

    float cap = (sinkInsts[i]->getInputCap() == 0)
                    ? options_->getSinkBufferInputCap()
                    : sinkInsts[i]->getInputCap();

    tarp.addSink(i, points[i].first, points[i].second, cap,
                 lp_target, hold_budget, nm);
  }

  // Hold budget summary (replaces per-sink CTS-0431 spam)
  logger_->info(CTS, 431,
      "TARP hold budget: {} direct, {} downstream-resolved, {} fail-closed (budget=0)",
      hb_direct, hb_downstream, hb_fail_closed);

  // Read timing graph CSV for setup slacks
  auto prof_csv_t0 = std::chrono::high_resolution_clock::now();
  const char* csv_env = std::getenv("CTS_FF_TIMING_GRAPH_CSV");
  if (csv_env) {
    std::ifstream csvFile(csv_env);
    if (csvFile.is_open()) {
      std::unordered_map<std::string, unsigned> nameToIdx;
      const auto& sinks = tarp.sinks();
      for (unsigned i = 0; i < sinks.size(); ++i)
        nameToIdx[sinks[i].name] = i;

      double T_period = 1.0;
      if (const char* e = std::getenv("CTS_CLOCK_PERIOD_NS"))
        T_period = std::stod(e);

      std::string line;
      std::getline(csvFile, line);  // skip header
      int edge_count = 0;
      while (std::getline(csvFile, line)) {
        std::vector<std::string> tokens;
        std::istringstream ls(line);
        std::string token;
        while (std::getline(ls, token, ',')) tokens.push_back(token);
        // CSV column layout (14 or 16 columns):
        //   0:from_ff, 1:to_ff, 2:slack_max, 3:slack_min, 4:arr_max, 5:arr_min,
        //   6:req_max, 7:req_min, 8:from_x, 9:from_y, 10:from_tier,
        //   11:to_x, 12:to_y, 13:to_tier [, 14:from_clock_net, 15:to_clock_net]
        // BUG FIX (V57_CLK2): from_clock_net/to_clock_net columns were appended later,
        // making the CSV 16 columns. The old "16-col legacy" branch (f_tier_idx=5,
        // t_tier_idx=8) incorrectly triggered, reading arrival_min_ns and from_x as tier
        // values — always mismatched — causing TARP to report 0 same-tier edges and pin
        // all FFs to Hilbert positions (timing-aware clustering completely disabled).
        // Fix: always use standard column positions (tier at col 10/13) for >=14 columns.
        if (tokens.size() < 14) continue;

        const int f_tier_idx = 10, t_tier_idx = 13;
        const int slack_max_idx = 2, slack_min_idx = 3;
        const int arr_max_idx = 4, arr_min_idx = 5;

        int f_tier = std::stoi(tokens[f_tier_idx]);
        int t_tier = std::stoi(tokens[t_tier_idx]);
        if (f_tier != t_tier) continue;  // skip cross-tier

        double slack_max = std::stod(tokens[slack_max_idx]);
        double slack_min = std::stod(tokens[slack_min_idx]);
        double arr_max = std::stod(tokens[arr_max_idx]);
        double arr_min = std::stod(tokens[arr_min_idx]);
        if (slack_max == 0 && slack_min == 0 && arr_max == 0 && arr_min == 0)
          continue;  // ghost edge

        auto it_a = nameToIdx.find(tokens[0]);
        auto it_b = nameToIdx.find(tokens[1]);
        if (it_a == nameToIdx.end() || it_b == nameToIdx.end()) continue;

        tarp.addTimingEdge(it_a->second, it_b->second, slack_max, T_period);
        ++edge_count;
      }
      logger_->info(CTS, 841, "TARP: {} same-tier edges from {}",
                    edge_count, csv_env);
    }
  }

  auto prof_csv_t1 = std::chrono::high_resolution_clock::now();
  double prof_csv_s = std::chrono::duration<double>(prof_csv_t1 - prof_csv_t0).count();
  logger_->info(CTS, 915, "Profile: TARP CSV parse {:.1f}s", prof_csv_s);

  auto prof_fm_t0 = std::chrono::high_resolution_clock::now();
  tarp.run();
  auto prof_fm_t1 = std::chrono::high_resolution_clock::now();
  double prof_fm_s = std::chrono::duration<double>(prof_fm_t1 - prof_fm_t0).count();
  logger_->info(CTS, 916, "Profile: TARP FM clustering {:.1f}s", prof_fm_s);

  // Convert clusters to leaf buffers
  std::vector<std::pair<float, float>> newSinkLocations;
  unsigned clusterCount = 0;

  for (const auto& cluster : tarp.clusters()) {
    if (cluster.empty()) continue;
    if (cluster.size() == 1) {
      unsigned orig = tarp.sinks()[cluster[0]].orig_idx;
      newSinkLocations.emplace_back(points[orig]);
      // Fix CTS-0492: register single-FF sink in mapLocationToSink_
      // (multi-FF clusters register at line 4725; single-FF was missing)
      const Point<double> sinkPt((double) points[orig].first,
                                 (double) points[orig].second);
      mapLocationToSink_[sinkPt] = const_cast<ClockInst*>(sinkInsts[orig]);
      ++clusterCount;
      continue;
    }

    std::vector<ClockInst*> clusterInsts;
    float xSum = 0, ySum = 0;
    for (unsigned si : cluster) {
      unsigned orig = tarp.sinks()[si].orig_idx;
      // Use sinkInsts directly — no map lookup needed, avoids float/double mismatch
      clusterInsts.push_back(const_cast<ClockInst*>(sinkInsts[orig]));
      xSum += points[orig].first;
      ySum += points[orig].second;
    }

    Point<double> center(xSum / cluster.size(), ySum / cluster.size());
    const int tier = cts3dDb_->getDominantTier(clusterInsts);
    const std::string buf = cts3dDb_->getBufferForTier(
        options_->getSinkBuffer(), tier);

    Point<double> loc = legalizeOneBuffer(center, buf);
    commitMoveLoc(center, loc);

    ClockInst& rootBuf = clock_.addClockBuffer(
        tierPrefix_ + "clkbuf_leaf_" + std::to_string(clusterCount), buf,
        loc.getX() * wireSegmentUnit_, loc.getY() * wireSegmentUnit_);
    cts3dDb_->setClockInstTier(rootBuf, tier);
    addFirstLevelSinkDriver(&rootBuf);

    ClockSubNet& net = clock_.addSubNet(
        tierPrefix_ + "clknet_leaf_" + std::to_string(clusterCount));
    net.addInst(rootBuf);
    for (ClockInst* ci : clusterInsts) net.addInst(*ci);
    net.setLeafLevel(true);

    const std::pair<float, float> bufPt(loc.getX(), loc.getY());
    newSinkLocations.emplace_back(bufPt);
    mapLocationToSink_[Point<double>(bufPt.first, bufPt.second)] = &rootBuf;

    if (cts3dDb_->hasSkewTargets()) {
      double sum = 0;
      for (unsigned si : cluster) sum += tarp.sinks()[si].lp_target;
      clusterBufTarget_[&rootBuf] = sum / cluster.size();
    }
    ++clusterCount;
  }

  topLevelSinksClustered_ = std::move(newSinkLocations);
  if (clusterCount) treeBufLevels_++;
  logger_->info(CTS, 844, "TARP: {} clusters, {} sinks",
                clusterCount, topLevelSinksClustered_.size());
}

// runTimingRefinement removed.
// FM min-cut with Hilbert initial partition replaces boundary swap.
// See TarpClustering.cpp for the FM implementation.

#if 0  // REMOVED: runTimingRefinement (replaced by TARP FM)
void HTreeBuilder::runTimingRefinement_REMOVED(
    std::vector<std::vector<unsigned>>& clusters,
    const std::vector<std::pair<float, float>>& points,
    const std::vector<const ClockInst*>& sinkInsts)
{
  const unsigned N = points.size();
  if (N == 0 || clusters.empty()) return;

  // Read params
  double slackThreshRatio = 0.3;
  if (const char* e = std::getenv("CTS_TARP_SLACK_THRESH"))
    slackThreshRatio = std::stod(e);

  // Get clock period (ns)
  double T_period = 1.0;
  if (const char* e = std::getenv("CTS_CLOCK_PERIOD_NS"))
    T_period = std::stod(e);
  const double slackThresh = slackThreshRatio * T_period;

  // Build FF name → point index map
  std::unordered_map<std::string, unsigned> nameToIdx;
  for (unsigned i = 0; i < N; ++i) {
    std::string nm = sinkInsts[i]->getName();
    const auto sp = nm.rfind('/');
    if (sp != std::string::npos) nm = nm.substr(0, sp);
    nameToIdx[nm] = i;
  }

  // Build point index → cluster index map
  std::vector<int> ffCluster(N, -1);
  for (unsigned c = 0; c < clusters.size(); ++c) {
    for (unsigned idx : clusters[c]) {
      if (idx < N) ffCluster[idx] = static_cast<int>(c);
    }
  }

  // Load critical timing edges from CSV
  const char* csv_env = std::getenv("PRE_CTS_FF_TIMING_GRAPH");
  if (!csv_env) {
    logger_->info(CTS, 875, "TARP refinement: no timing graph CSV, skipping");
    return;
  }

  struct CritEdge {
    unsigned ff_a, ff_b;
    double slack;
  };
  std::vector<CritEdge> critEdges;

  std::ifstream csvFile(csv_env);
  if (!csvFile.is_open()) {
    logger_->warn(CTS, 876, "TARP refinement: cannot open {}", csv_env);
    return;
  }

  std::string line;
  std::getline(csvFile, line);  // skip header
  while (std::getline(csvFile, line)) {
    std::stringstream ss(line);
    std::string token;
    std::vector<std::string> tokens;
    while (std::getline(ss, token, ',')) tokens.push_back(token);
    if (tokens.size() < 14) continue;

    const std::string& from_ff = tokens[0];
    const std::string& to_ff = tokens[1];
    double slack_max = std::stod(tokens[2]);
    int from_tier = std::stoi(tokens[10]);
    int to_tier = std::stoi(tokens[13]);

    // Same-tier only
    if (from_tier != to_tier) continue;
    // Ghost edge filter
    double slack_min = std::stod(tokens[3]);
    double arr_max = std::stod(tokens[4]);
    double arr_min = std::stod(tokens[5]);
    if (slack_max == 0 && slack_min == 0 && arr_max == 0 && arr_min == 0) continue;
    // Slack threshold filter
    if (slack_max > slackThresh) continue;

    auto it_a = nameToIdx.find(from_ff);
    auto it_b = nameToIdx.find(to_ff);
    if (it_a == nameToIdx.end() || it_b == nameToIdx.end()) continue;

    critEdges.push_back({it_a->second, it_b->second, slack_max});
  }

  // Count cross-cluster critical pairs
  unsigned crossBefore = 0;
  for (const auto& e : critEdges) {
    if (ffCluster[e.ff_a] != ffCluster[e.ff_b]) crossBefore++;
  }

  logger_->info(CTS, 877,
      "TARP refinement: {} critical edges (slack<{:.3f}ns), {} cross-cluster ({:.1f}%)",
      critEdges.size(), slackThresh,
      crossBefore, critEdges.empty() ? 0.0 : 100.0 * crossBefore / critEdges.size());

  if (crossBefore == 0) {
    logger_->info(CTS, 878, "TARP refinement: no cross-cluster critical pairs, skipping");
    return;
  }

  // Compute max diameter from best clustering (approximate from existing clusters)
  double maxDiameter = 0;
  for (const auto& cl : clusters) {
    if (cl.size() <= 1) continue;
    float xn = points[cl[0]].first, xx = xn;
    float yn = points[cl[0]].second, yx = yn;
    for (unsigned idx : cl) {
      xn = std::min(xn, points[idx].first); xx = std::max(xx, points[idx].first);
      yn = std::min(yn, points[idx].second); yx = std::max(yx, points[idx].second);
    }
    double diam = (xx - xn) + (yx - yn);
    maxDiameter = std::max(maxDiameter, diam);
  }
  // Allow 20% diameter expansion for refinement
  maxDiameter *= 1.2;

  // Build per-FF critical edge count (how many critical edges this FF participates in)
  std::vector<unsigned> ffCritCount(N, 0);
  for (const auto& e : critEdges) {
    if (ffCluster[e.ff_a] != ffCluster[e.ff_b]) {
      ffCritCount[e.ff_a]++;
      ffCritCount[e.ff_b]++;
    }
  }

  // Greedy swap: for each cross-cluster critical edge, try moving the FF
  // with fewer critical edges to the other's cluster
  unsigned swaps = 0;
  const unsigned maxSwaps = N / 4;  // cap at 25% of FFs

  // Sort critical edges by slack (tightest first — most important to fix)
  auto sortedEdges = critEdges;
  std::sort(sortedEdges.begin(), sortedEdges.end(),
            [](const CritEdge& a, const CritEdge& b) { return a.slack < b.slack; });

  for (const auto& e : sortedEdges) {
    if (swaps >= maxSwaps) break;

    int ca = ffCluster[e.ff_a];
    int cb = ffCluster[e.ff_b];
    if (ca == cb) continue;  // already co-clustered (maybe by previous swap)

    // Pick the FF to move: the one with fewer critical edges in its current cluster
    // (less disruptive to move)
    unsigned mover, target_cluster;
    if (ffCritCount[e.ff_a] <= ffCritCount[e.ff_b]) {
      mover = e.ff_a;
      target_cluster = cb;
    } else {
      mover = e.ff_b;
      target_cluster = ca;
    }

    int src_cluster = ffCluster[mover];

    // Don't empty a cluster
    if (clusters[src_cluster].size() <= 2) continue;

    // Check diameter constraint after swap
    float xn = points[mover].first, xx = xn;
    float yn = points[mover].second, yx = yn;
    for (unsigned idx : clusters[target_cluster]) {
      xn = std::min(xn, points[idx].first); xx = std::max(xx, points[idx].first);
      yn = std::min(yn, points[idx].second); yx = std::max(yx, points[idx].second);
    }
    double newDiam = (xx - xn) + (yx - yn);
    if (newDiam > maxDiameter) continue;  // would violate spatial constraint

    // Execute swap
    // Remove from source cluster
    auto& srcCluster = clusters[src_cluster];
    srcCluster.erase(std::remove(srcCluster.begin(), srcCluster.end(), mover),
                     srcCluster.end());
    // Add to target cluster
    clusters[target_cluster].push_back(mover);
    ffCluster[mover] = target_cluster;
    swaps++;
  }

  // Count cross-cluster after refinement
  unsigned crossAfter = 0;
  for (const auto& e : critEdges) {
    if (ffCluster[e.ff_a] != ffCluster[e.ff_b]) crossAfter++;
  }

  logger_->info(CTS, 879,
      "TARP refinement: {} swaps, cross-cluster {} -> {} ({:.1f}% reduction)",
      swaps, crossBefore, crossAfter,
      crossBefore > 0 ? 100.0 * (crossBefore - crossAfter) / crossBefore : 0.0);
}
#endif  // REMOVED: runTimingRefinement

}  // namespace cts
