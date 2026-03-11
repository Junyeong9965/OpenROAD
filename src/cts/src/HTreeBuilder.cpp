// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "HTreeBuilder.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "Clustering.h"
#include "Cts3DDatabase.h"  // JYJ (2026-02-06) Added for 3D tier management
#include "SinkClustering.h"
#include "TechChar.h"
#include "TreeBuilder.h"
#include "Util.h"
#include "odb/db.h"
#include "odb/isotropy.h"
#include "utl/Logger.h"

namespace cts {

using utl::CTS;

// JYJ (2026-02-06) Removed anonymous namespace containing hasSuffix() and
// mapBufferMasterToTier() — moved to Cts3DDatabase::hasSuffix() and
// Cts3DDatabase::getBufferForTier()

Point<double> HTreeBuilder::legalizeOneBuffer(Point<double> bufferLoc,
                                              const std::string& bufferName)
{
  Point<double> legalLoc
      = TreeBuilder::legalizeOneBuffer(bufferLoc, bufferName);
  return resolveLocationCollision(legalLoc);
}

// JYJ (2026-02-21) Step 3: Skew-aware clustering via skew target penalty
// HB penalty removed — vias (including HB) don't get distance penalty in 2D CTS
// either, and STA already accounts for HB RC in timing. Clustering should remain
// purely spatial + skew-target-driven.
double HTreeBuilder::computeDist(const Point<double>& x, const Point<double>& y)
{
  double baseDist = TreeBuilder::computeDist(x, y);

  // Skew target penalty: group FFs with similar arrival targets
  if (skewTargetBeta_ > 0.0 && cts3dDb_ != nullptr
      && cts3dDb_->hasSkewTargets()) {
    auto itX = mapLocationToSink_.find(x);
    auto itY = mapLocationToSink_.find(y);
    if (itX != mapLocationToSink_.end() && itY != mapLocationToSink_.end()
        && itX->second != nullptr && itY->second != nullptr) {
      // V32c: use getClockInstTarget() which handles both cluster buffers and FFs
      const double targetX = getClockInstTarget(itX->second);
      const double targetY = getClockInstTarget(itY->second);
      baseDist += skewTargetBeta_ * std::abs(targetX - targetY);
    }
  }

  return baseDist;
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
      const Point<double> normLocation((float) inst.getX() / wireSegmentUnit_,
                                       (float) inst.getY() / wireSegmentUnit_);
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

  // JYJ (2026-02-09) Tier-aware clustering: detect cross-tier clusters
  // Calculate HB equivalent distance for cross-tier penalty awareness
  const double wire_res_per_unit = cts3dDb_->getResPerDBU(0) * wireSegmentUnit_;
  const double wire_cap_per_unit = cts3dDb_->getCapPerDBU(0) * wireSegmentUnit_;
  const double hb_equivalent_dist = cts3dDb_->getHbtEquivalentDistance(
      wire_res_per_unit, wire_cap_per_unit);

  // Note: Current k-means clustering uses geometric distance only.
  // TODO Phase 3: Modify SinkClustering/CKMeans to inject hb_equivalent_dist
  // as penalty when computing distance between sinks on different tiers.
  // This would actively discourage cross-tier clusters during k-means.

  std::vector<std::pair<float, float>> newSinkLocations;

  // JYJ (2026-02-25) V37: Helper lambda to create one leaf buffer for a group of FFs.
  // This is factored out so that sub-clustering can call it multiple times per cluster.
  int splitCount_2way = 0, splitCount_3way = 0;  // V37 stats
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

    const char* baseName = secondLevel ? "clkbuf_leaf2_" : "clkbuf_leaf_";
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

    const char* netBaseName = secondLevel ? "clknet_leaf2_" : "clknet_leaf_";
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

    // V32c: store mean LP skew target for this leaf buffer
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
          logger_->error(CTS, 79, "Sink not found.");
        }
        clusterClockInsts.push_back(mapLocationToSink_[mapPoint]);
      }

      // JYJ (2026-02-09) Detect cross-tier clustering
      int tier0_count = 0, tier1_count = 0;
      for (const ClockInst* inst : clusterClockInsts) {
        const int inst_tier = cts3dDb_->getInstTier(inst->getDbInst());
        if (inst_tier == 0) tier0_count++;
        else if (inst_tier == 1) tier1_count++;
      }
      if (tier0_count > 0 && tier1_count > 0) {
        logger_->warn(CTS, 364,
                      "Cross-tier cluster detected: cluster={}, tier0={}, tier1={}, "
                      "dominant_tier={}, HB_penalty={:.1f}um",
                      clusterCount, tier0_count, tier1_count,
                      cts3dDb_->getDominantTier(clusterClockInsts),
                      hb_equivalent_dist);
      }

      // JYJ (2026-03-03) V48: Lambda for V37 target-aware split + leaf buffer
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
                              "V37 sub-cluster {}.{}: {} sinks, "
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

      // JYJ (2026-03-01) V48: Tier-aware cluster split.
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
            "V48 tier split cluster {}: bottom={}, upper={}",
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

  // V37: Log sub-clustering statistics
  if (splitCount_2way > 0 || splitCount_3way > 0) {
    logger_->info(CTS, 414,
                  "V37 target-aware sub-clustering: {} 2-way splits, "
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

  // V37-debug: Verify every entry in topLevelSinksClustered_ exists in mapLocationToSink_
  int missingCount = 0;
  for (size_t si = 0; si < topLevelSinksClustered_.size(); ++si) {
    const auto& sinkPair = topLevelSinksClustered_[si];
    Point<double> key(sinkPair.first, sinkPair.second);
    if (mapLocationToSink_.find(key) == mapLocationToSink_.end()) {
      logger_->warn(CTS, 418,
                    "V37-debug: topLevelSinksClustered_[{}] = ({:.10f}, {:.10f}) "
                    "NOT in mapLocationToSink_ (size={})",
                    si, sinkPair.first, sinkPair.second,
                    mapLocationToSink_.size());
      ++missingCount;
    }
  }
  if (missingCount > 0) {
    logger_->warn(CTS, 419,
                  "V37-debug: {} / {} sinks missing from mapLocationToSink_!",
                  missingCount, topLevelSinksClustered_.size());
  } else {
    logger_->info(CTS, 427,
                  "V37-debug: All {} sinks verified in mapLocationToSink_ (mapSize={})",
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
    // Simulate the float conversion that happens when storing in
    // newSinkLocations. This is needed because the existing logic uses the
    // double->float conversion for the center location.
    Point<double> checkKey((float) resolvedLocation.getX(),
                           (float) resolvedLocation.getY());
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

// JYJ (2026-02-06) Removed getDominantTierFromInsts, getDominantTierFromSinkLocs,
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

  // JYJ (2026-02-21) Step 3: Initialize skew-aware clustering weight from env
  skewTargetBeta_ = 0.0;
  if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()) {
    const char* beta_env = std::getenv("CTS_SKEW_TARGET_BETA");
    if (beta_env != nullptr) {
      skewTargetBeta_ = std::atof(beta_env);
    } else {
      skewTargetBeta_ = 5000.0;  // default
    }
    logger_->info(CTS, 400,
                  "Skew-aware clustering enabled: beta={:.1f}, {} targets loaded",
                  skewTargetBeta_, cts3dDb_->getSkewTargetCount());
  }

  // JYJ (2026-02-23) V31: Initialize useful-skew wire adjustment parameters
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
      wireSkewScale_ = 1.0;  // default: enabled when skew targets are loaded
    }
    if (wdpu_env != nullptr) {
      wireDelayPerUnit_ = std::atof(wdpu_env) * 0.001;  // ps/um → ns/um
    }
    if (wireSkewScale_ > 0.0) {
      logger_->info(CTS, 405,
                    "V31 useful-skew wire adjustment enabled: scale={:.2f}, "
                    "wireDelay={:.3f}ps/um",
                    wireSkewScale_, wireDelayPerUnit_ * 1000.0);
    }
  }

  // JYJ (2026-02-27) V42-fix: Read target-split env vars BEFORE initSinkRegion(),
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

  // JYJ (2026-03-01) V48: Tier-aware cluster split
  if (const char* e = std::getenv("CTS_ENABLE_TIER_SPLIT"))
    enableTierSplit_ = (std::atoi(e) != 0);
  logger_->info(CTS, 376, "V48 tier split: {}", enableTierSplit_ ? "ON" : "OFF");

  initSinkRegion();

  for (int level = 1; level <= clockTreeMaxDepth_; ++level) {
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

  // JYJ (2026-02-21) Step 4: Compute per-branch delay targets before tree construction
  // JYJ (2026-02-23) V32b: also read N-buffer chain parameters (singleBufDelay_, maxLeafDelayBufs_)
  if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()) {
    const char* thr_env = std::getenv("CTS_DELAY_TARGET_THRESHOLD");
    delayTargetThreshold_ = (thr_env != nullptr) ? std::atof(thr_env) : 0.010;

    // V32b: single delay-buffer delay (ns) and max chain length per leaf
    if (const char* e = std::getenv("CTS_LEAF_BUF_DELAY_NS"))
      singleBufDelay_ = std::atof(e);
    if (const char* e = std::getenv("CTS_MAX_LEAF_DELAY_BUFS"))
      maxLeafDelayBufs_ = std::atoi(e);

    // JYJ (2026-02-25) V37: Target-aware sub-clustering parameters
    // NOTE: env var reading moved to before initSinkRegion() (V42-fix).
    // Only logging remains here.
    if (enableTargetSplit_) {
      logger_->info(CTS, 410,
                    "V37 target-aware sub-clustering enabled: "
                    "2way_threshold={:.1f}ps, 3way_threshold={:.1f}ps",
                    splitThreshold2wayNs_ * 1000.0,
                    splitThreshold3wayNs_ * 1000.0);
    }

    // JYJ (2026-02-26) V39: Per-FF relay buffer parameters
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
      if (const char* csv = std::getenv("CTS_TIMING_GRAPH_CSV"))
        loadPerFfHoldBudgets(csv);
      // JYJ V40: Also load IO edge hold budgets (PI→FF hold slack)
      if (const char* ioCsv = std::getenv("CTS_IO_TIMING_CSV"))
        loadIoHoldBudgets(ioCsv);
      // JYJ (2026-02-28) V47: Elmore wire RC parameters for x_useful positioning
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
                    "V47 per-FF relay: bufDelay={:.3f}ns maxRelay={} "
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

    // JYJ (2026-03-05) V49: Grouped Delay Chain parameters.
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
        if (const char* csv = std::getenv("CTS_TIMING_GRAPH_CSV"))
          loadPerFfHoldBudgets(csv);
        if (const char* ioCsv = std::getenv("CTS_IO_TIMING_CSV"))
          loadIoHoldBudgets(ioCsv);
      }
      // V50: Cluster-Uniform Depth mode — one chain per cluster, depth = cluster median.
      enableClusterUniform_ = false;
      if (const char* e = std::getenv("CTS_GROUPED_DELAY_CLUSTER_UNIFORM"))
        enableClusterUniform_ = (std::atoi(e) != 0);
      logger_->info(CTS, 490,
                    "V50 grouped delay: maxDepth={} minDelta={:.3f}ns "
                    "bufDelay={:.3f}ns holdGuard={:.3f}ns holdBudgets={} "
                    "clusterUniform={}",
                    groupedDelayMaxDepth_, groupedDelayMinDelta_,
                    perFfBufDelay_, perFfHoldGuard_,
                    static_cast<int>(perFfHoldBudget_.size()),
                    enableClusterUniform_ ? 1 : 0);
    } else if (enableGroupedDelay_) {
      logger_->warn(CTS, 491,
                    "CTS_ENABLE_GROUPED_DELAY=1 but no skew targets; "
                    "falling back to balanced CTS");
      enableGroupedDelay_ = false;
    }

    // V39: Intermediate delay buffer gate (V37 default on, V39 default off)
    enableMidDelayBufs_ = false;
    if (const char* e = std::getenv("CTS_ENABLE_MID_DELAY_BUFS"))
      enableMidDelayBufs_ = (std::atoi(e) != 0);

    // JYJ (2026-02-26) V40: CKMeans target-aware branching weight.
    // When > 0, H-tree 2-way split (CKMeans) groups FFs with similar targets.
    ckmeansTargetBeta_ = 0.0;
    if (const char* e = std::getenv("CTS_CKMEANS_TARGET_BETA"))
      ckmeansTargetBeta_ = std::atof(e);
    if (ckmeansTargetBeta_ > 0.0) {
      logger_->info(CTS, 428,
          "V40 CKMeans target-aware branching: beta={:.1f}",
          ckmeansTargetBeta_);
    }

    computeBranchDelayTargets();
  }

  createClockSubNets();
  // clang-format off
  debugPrint(logger_, CTS, "legalizer", 3, "Htree file {} has been generated",
             plotHTree());
  debugPrint(logger_, CTS, "legalizer", 3, "Run 'obsAwareCts.py cts.clk.buffer'"
	     "to produce cts.clk.buffer.png");
  // clang-format on
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

void HTreeBuilder::initTopLevelSinks(
    std::vector<std::pair<float, float>>& sinkLocations,
    std::vector<const ClockInst*>& sinkInsts)
{
  sinkLocations.clear();
  clock_.forEachSink([&](const ClockInst& sink) {
    sinkLocations.emplace_back((float) sink.getX() / wireSegmentUnit_,
                               (float) sink.getY() / wireSegmentUnit_);
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

  // JYJ (2026-02-26) V40: Pass per-sink LP targets to CKMeans.
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

  // JYJ (2026-02-23) V31: Useful-skew wire length adjustment (Fishburn 1990)
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
          // V32c: unified target accessor (cluster buffer or individual FF)
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
      // JYJ V41-fix BUG#6: Widened clamp from 50% to 80% of branch length.
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
  // End JYJ V31

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

// JYJ (2026-02-23) V32c: unified LP skew target accessor for any ClockInst.
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

// JYJ (2026-02-21) Step 4: Compute mean arrival target per branch.
// JYJ (2026-02-25) V37: Extended to ALL levels (not just leaf).
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

  // V32c-debug: per-branch LP target statistics
  struct BranchDetail {
    int level;
    unsigned idx;
    int count;
    double mean, bmin, bmax, stddev;
  };
  std::vector<BranchDetail> branchDetails;

  // V37: Compute per-branch mean target at ALL levels (not just leaf).
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
        auto it = mapLocationToSink_.find(loc);
        if (it == mapLocationToSink_.end() || it->second == nullptr) {
          continue;
        }
        const double tgt = getClockInstTarget(it->second);
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

  // V32c-debug CTS-407: per-branch LP target distribution detail (leaf level only)
  for (const auto& d : branchDetails) {
    if (d.level != leafLevelIdx) continue;
    const double delta = d.mean - globalMeanTarget_;
    int plannedN = 0;
    if (delta > 0.0 && singleBufDelay_ > 1e-9) {
      plannedN = std::max(0, std::min(
          static_cast<int>(std::round(delta / singleBufDelay_)),
          maxLeafDelayBufs_));
    }
    logger_->info(CTS, 407,
                  "  Leaf branch {:2d}: sinks={}, "
                  "LP=[min={:.4f} mean={:.4f} max={:.4f} spread={:.4f} "
                  "stddev={:.4f}]ns, delta={:.4f}ns -> planned_N={}",
                  d.idx, d.count,
                  d.bmin, d.mean, d.bmax, d.bmax - d.bmin,
                  d.stddev, delta, plannedN);
  }
}

// JYJ (2026-02-26) V39: Load per-FF hold budgets from timing graph CSV.
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
                "V39 hold budgets loaded: {} capture FFs from {} edges",
                perFfHoldBudget_.size(), edgeCount);
}

// JYJ (2026-02-26) V40: Load IO edge hold budgets from IO timing CSV.
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
                "V40 IO hold budgets: {} PI→FF edges, {} new FFs, "
                "{} tightened (total hold budget FFs: {})",
                ioEdgeCount, newCount, mergedCount,
                perFfHoldBudget_.size());
}

void HTreeBuilder::createClockSubNets()
{
  Point<double> center = sinkRegion_.getCenter();
  // JYJ (2026-02-06) Replaced with Cts3DDatabase calls
  const int root_tier = cts3dDb_->getDominantTierFromClock(clock_);
  const std::string root_buffer
      = cts3dDb_->getBufferForTier(options_->getRootBuffer(), root_tier);
  logger_->info(utl::CTS, 315,
                "3D-CTS root buffer: tier={}, buffer={}",
                root_tier, root_buffer);
  Point<double> legalCenter = legalizeOneBuffer(center, root_buffer);
  sinkRegion_.setCenter(legalCenter);
  commitMoveLoc(center, legalCenter);
  const int centerX = legalCenter.getX() * wireSegmentUnit_;
  const int centerY = legalCenter.getY() * wireSegmentUnit_;

  ClockInst& rootBuffer
      = clock_.addClockBuffer("clkbuf_0", root_buffer, centerX, centerY);

  if (topBufferName_.empty()) {
    topBufferName_ = rootBuffer.getName();
  }

  // clang-format off
  if (center != legalCenter) {
    debugPrint(logger_, CTS, "legalizer", 2, "createClockSubNets: "
	       "root clkbuf_0: {} => {}", center, legalCenter);
  } else {
    debugPrint(logger_, CTS, "legalizer", 2, "createClockSubNets: "
	       "root clkbuf_0: {}", center);
  }
  // clang-format on

  addTreeLevelBuffer(&rootBuffer);
  ClockSubNet& rootClockSubNet = clock_.addSubNet("clknet_0");
  rootClockSubNet.setTreeLevel(0);
  rootClockSubNet.addInst(rootBuffer);
  treeBufLevels_++;

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
    // JYJ (2026-02-06) Replaced with Cts3DDatabase calls
    const int branch_tier = cts3dDb_->getDominantTier(
        topLevelTopology.getBranchSinksLocations(idx), mapLocationToSink_);
    const std::string branch_root_buffer = cts3dDb_->getBufferForTier(
        options_->getRootBuffer(), branch_tier);
    logger_->info(utl::CTS, 316,
                  "3D-CTS branch L1 idx={}: tier={}, buffer={}",
                  idx, branch_tier, branch_root_buffer);

    // JYJ (2026-02-09, updated 2026-03-10) Cross-tier HB delay consideration
    // Elmore 50% delay: t_via = 0.693 * R_HB * C_HB (negligible: ~0 with C=0)
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

    SegmentBuilder builder("clkbuf_1_" + std::to_string(idx) + "_",
                           "clknet_1_" + std::to_string(idx) + "_",
                           legalCenter,  // center may have moved, don't use
                                         // sinkRegion_.getCenter()
                           legalBranchPoint,
                           topLevelTopology.getWireSegments(),
                           clock_,
                           rootClockSubNet,
                           *techChar_,
                           wireSegmentUnit_,
                           this,
                           db_,
                           branch_tier);
    if (!options_->getTreeBuffer().empty()) {
      // JYJ (2026-02-06) Replaced mapBufferMasterToTier with Cts3DDatabase
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

    // JYJ (2026-02-25) V37: Insert intermediate delay buffers at Level 1.
    // If this branch's subtree mean target > level globalMean, add N delay
    // buffers in series to shift arrival for the entire subtree.
    ClockSubNet* branchDrivingSub = builder.getDrivingSubNet();
    // JYJ (2026-02-26) V39: gate with enableMidDelayBufs_ (default off in V39)
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
            const std::string bufName = "clkbuf_mid_dly_1_"
                + std::to_string(idx) + "_" + std::to_string(b);
            ClockInst& buf = clock_.addClockBuffer(
                bufName, delayBufMaster, bx, by);
            cts3dDb_->setClockInstTier(buf, branch_tier);
            addTreeLevelBuffer(&buf);
            curSub->addInst(buf);
            ClockSubNet& nextSub = clock_.addSubNet(
                "clknet_mid_dly_1_" + std::to_string(idx)
                + "_" + std::to_string(b));
            nextSub.addInst(buf);
            curSub = &nextSub;
          }
          branchDrivingSub = curSub;
          logger_->info(CTS, 411,
                        "V37 intermediate delay L1 branch {}: "
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

      // JYJ (2026-02-06) Replaced with Cts3DDatabase calls
      const int branch_tier = cts3dDb_->getDominantTier(
          topology.getBranchSinksLocations(idx), mapLocationToSink_);
      const std::string branch_root_buffer = cts3dDb_->getBufferForTier(
          options_->getRootBuffer(), branch_tier);
      logger_->info(utl::CTS, 317,
                    "3D-CTS branch L{} idx={}: tier={}, buffer={}",
                    levelIdx+1, idx, branch_tier, branch_root_buffer);

      // JYJ (2026-02-09, updated 2026-03-10) Cross-tier HB delay consideration
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
		   "{} : {} => {}", levelIdx+1, "clkbuf_" +
		   std::to_string(levelIdx+1) + "_" + std::to_string(idx) + "_",
		   branchPoint, legalBranchPoint);
      } else {
	debugPrint(logger_, CTS, "legalizer", 2, "createClockSubNets level {} "
		   "{} : {}", levelIdx+1, "clkbuf_" +
		   std::to_string(levelIdx+1) + "_" + std::to_string(idx) + "_",
		   branchPoint);
      }
      // clang-format on

      SegmentBuilder builder("clkbuf_" + std::to_string(levelIdx + 1) + "_"
                                 + std::to_string(idx) + "_",
                             "clknet_" + std::to_string(levelIdx + 1) + "_"
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
        // JYJ (2026-02-06) Replaced mapBufferMasterToTier with Cts3DDatabase
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

      // JYJ (2026-02-25) V37: Insert intermediate delay buffers at Level 2+.
      // Skip the leaf level (handled by existing N-chain code below).
      const int leafLevel
          = static_cast<int>(topologyForEachLevel_.size()) - 1;
      ClockSubNet* branchDrivingSub = builder.getDrivingSubNet();
      // JYJ (2026-02-26) V39: gate with enableMidDelayBufs_ (default off in V39)
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
              const std::string bufName = "clkbuf_mid_dly_"
                  + std::to_string(levelIdx + 1) + "_"
                  + std::to_string(idx) + "_" + std::to_string(b);
              ClockInst& buf = clock_.addClockBuffer(
                  bufName, delayBufMaster, bx, by);
              cts3dDb_->setClockInstTier(buf, branch_tier);
              addTreeLevelBuffer(&buf);
              curSub->addInst(buf);
              ClockSubNet& nextSub = clock_.addSubNet(
                  "clknet_mid_dly_" + std::to_string(levelIdx + 1)
                  + "_" + std::to_string(idx) + "_" + std::to_string(b));
              nextSub.addInst(buf);
              curSub = &nextSub;
            }
            branchDrivingSub = curSub;
            logger_->info(CTS, 412,
                          "V37 intermediate delay L{} branch {}: "
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

  // V32c-debug: LP-CTS gap tracking accumulators (reported in CTS-408 after loop).
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

        // JYJ (2026-02-21) Step 4: Insert delay buffers if branch LP target > globalMean
        // JYJ (2026-02-23) V32b: N-buffer chain based on ABSOLUTE LP-TNS target.
        //   Problem (found after V32b run): Absolute targeting causes hold violations.
        //   Branch differentials up to 36-48ps (N=1 vs N=4) violated hold on short paths.
        //   LP guaranteed hold only for per-FF delays; branch-mean averaging breaks this.
        // JYJ (2026-02-23) V32c-fix: DELTA-BASED N-buffer chain (relative to globalMean).
        //   N = round(max(0, branchMean - globalMean) / singleBufDelay_)
        //   Limits max inter-branch differential to (maxBranchDelta - 0) ≈ 12-16ps (1-2 bufs).
        //   LP-TNS result (asap7/aes): globalMean=26.6ps, max delta=15.8ps → N≤1 per branch.
        //   Chain: subNet → buf_0 → sub_0 → ... → FFs
        // V32c-debug: branchCnt hoisted to lambda scope for CTS-408 gap tracking
        double absTgt = 0.0;
        int branchCnt = 0;
        if (cts3dDb_ != nullptr && cts3dDb_->hasSkewTargets()
            && singleBufDelay_ > 1e-9) {
          const auto& bSinkLocs = leafTopology.getBranchSinksLocations(idx);
          double branchSum = 0.0;
          for (const auto& bLoc : bSinkLocs) {
            auto bIt = mapLocationToSink_.find(bLoc);
            if (bIt == mapLocationToSink_.end() || bIt->second == nullptr) {
              continue;
            }
            // V32c: unified target accessor (cluster buffer or individual FF)
            branchSum += getClockInstTarget(bIt->second);
            ++branchCnt;
          }
          if (branchCnt > 0) {
            absTgt = branchSum / branchCnt;
          }
        }

        // JYJ (2026-03-05) V49: Grouped Delay Chain.
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
          // V49-fix2: static counter for globally unique buffer names.
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
          // V50: When enableClusterUniform_=true, all FFs in the cluster share
          // a single depth k_c computed from the cluster median LP target delta.
          // This creates exactly 1 chain per cluster (vs up to MAX_DEPTH chains
          // in V49), reducing buffer count by ~60% with comparable timing.
          // Hold safety: k_c is capped by the minimum hold budget across all FFs.
          auto collectFFs = [&](ClockInst* driverInst, ClockSubNet* driverSubNet,
                                std::vector<FfGroupInfo>& ffInfos,
                                std::vector<ClockInst*>& memberFFs) {
            if (enableClusterUniform_) {
              // V50: Cluster-Uniform mode — compute median delta, assign to all FFs.
              // Pass 1: collect deltas and hold budgets
              std::vector<double> deltas;
              deltas.reserve(memberFFs.size());
              for (ClockInst* ff : memberFFs) {
                std::string ffName = getFFName(ff);
                double ffDelta = cts3dDb_->getSkewTarget(ffName) - globalMeanTarget_;
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

              // Hold safety: cap by minimum hold budget across all FFs in cluster
              if (clusterDepth > 0) {
                for (ClockInst* ff : memberFFs) {
                  std::string ffName = getFFName(ff);
                  auto hIt = perFfHoldBudget_.find(ffName);
                  if (hIt != perFfHoldBudget_.end()) {
                    while (clusterDepth > 0
                        && hIt->second - clusterDepth * perFfBufDelay_
                               < perFfHoldGuard_) {
                      --clusterDepth;
                    }
                  }
                }
                if (clusterDepth == 0) ++branchHoldSkip;
              }

              // Pass 2: assign uniform depth to all FFs
              for (size_t i = 0; i < memberFFs.size(); ++i) {
                ffInfos.push_back({memberFFs[i], deltas[i], clusterDepth});
              }
            } else {
              // V49: Per-FF depth (original behavior)
              for (ClockInst* ff : memberFFs) {
                std::string ffName = getFFName(ff);
                double ffTarget = cts3dDb_->getSkewTarget(ffName);
                double ffDelta = ffTarget - globalMeanTarget_;

                // Quantize to buffer delay steps
                int depth = 0;
                if (ffDelta > groupedDelayMinDelta_) {
                  depth = static_cast<int>(
                      std::round(ffDelta / perFfBufDelay_));
                  depth = std::max(0, std::min(depth, groupedDelayMaxDepth_));
                }

                // Hold safety: cap depth by hold budget
                if (depth > 0) {
                  auto hIt = perFfHoldBudget_.find(ffName);
                  if (hIt != perFfHoldBudget_.end()) {
                    while (depth > 0
                        && hIt->second - depth * perFfBufDelay_ < perFfHoldGuard_) {
                      --depth;
                    }
                    if (depth == 0) ++branchHoldSkip;
                  }
                }

                ffInfos.push_back({ff, ffDelta, depth});
              }
            }
          };

          // Build grouped delay chains for collected FFs
          // V49-fix3: No local numSinks here. Uses outer numSinks (line 3002)
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

                // V49a: Filter out cross-tier minority FFs before building the chain.
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
                ClockSubNet* cur = driverSubNet;
                for (int d = 0; d < depth; ++d) {
                  double frac = static_cast<double>(d + 1) / (depth + 1);
                  int bx = driverInst->getX()
                      + static_cast<int>((cx - driverInst->getX()) * frac);
                  int by = driverInst->getY()
                      + static_cast<int>((cy - driverInst->getY()) * frac);
                  // V50-Legal (JYJ 2026-03-10): Grid-snap grouped delay buffer via
                  // legalizeOneBuffer(), same as all H-tree standard buffers.
                  // Raw centroid interpolation produces off-grid DBU positions → DPL-0036.
                  // Conversion: DBU → techChar units → legalizeOneBuffer → DBU.
                  {
                    Point<double> rawLoc(static_cast<double>(bx) / wireSegmentUnit_,
                                        static_cast<double>(by) / wireSegmentUnit_);
                    Point<double> legalLoc = legalizeOneBuffer(rawLoc, grpBufMaster);
                    bx = static_cast<int>(legalLoc.getX() * wireSegmentUnit_);
                    by = static_cast<int>(legalLoc.getY() * wireSegmentUnit_);
                  }

                  // V49-fix2: use myGrpId (global counter) instead of idx
                  // to avoid collision across sub-trees.
                  const std::string bName = "clkbuf_grp_"
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
                      "clknet_grp_" + std::to_string(myGrpId) + "_c"
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

          // Iterate sinks in this branch (same pattern as per-FF relay)
          const std::vector<Point<double>>& sinkLocs
              = leafTopology.getBranchSinksLocations(idx);
          for (const auto& loc : sinkLocs) {
            // Epsilon-tolerant sink lookup (V37 CTS-0080 fix)
            auto sinkIt = mapLocationToSink_.find(loc);
            if (sinkIt == mapLocationToSink_.end()) {
              constexpr double eps = 1e-6;
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
              logger_->error(CTS, 492,
                  "V49 grouped: Sink not found at ({:.10f}, {:.10f})",
                  loc.getX(), loc.getY());
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
                buildGroupChains(sinkInst, leafSubNet, ffInfos);
              } else {
                ++branchDirectFFs;
              }
            } else {
              // Individual FF (no clustering): treat as single-FF group
              std::string ffName = getFFName(sinkInst);
              double ffTarget = cts3dDb_->getSkewTarget(ffName);
              double ffDelta = ffTarget - globalMeanTarget_;
              int depth = 0;
              if (ffDelta > groupedDelayMinDelta_) {
                depth = static_cast<int>(
                    std::round(ffDelta / perFfBufDelay_));
                depth = std::max(0, std::min(depth, groupedDelayMaxDepth_));
              }
              if (depth > 0) {
                auto hIt = perFfHoldBudget_.find(ffName);
                if (hIt != perFfHoldBudget_.end()) {
                  while (depth > 0
                      && hIt->second - depth * perFfBufDelay_ < perFfHoldGuard_) {
                    --depth;
                  }
                  if (depth == 0) ++branchHoldSkip;
                }
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
                ClockSubNet* cur = subNet;
                for (int d = 0; d < depth; ++d) {
                  double frac = static_cast<double>(d + 1) / (depth + 1);
                  int bx = cur->getDriver()->getX()
                      + static_cast<int>(
                          (sinkInst->getX() - cur->getDriver()->getX()) * frac);
                  int by = cur->getDriver()->getY()
                      + static_cast<int>(
                          (sinkInst->getY() - cur->getDriver()->getY()) * frac);
                  // V50-Legal (JYJ 2026-03-10): Grid-snap per-FF grouped delay buffer
                  // via legalizeOneBuffer(), same as all H-tree standard buffers.
                  // Raw interpolation produces off-grid DBU positions → DPL-0036.
                  // Conversion: DBU → techChar units → legalizeOneBuffer → DBU.
                  {
                    Point<double> rawLoc(static_cast<double>(bx) / wireSegmentUnit_,
                                        static_cast<double>(by) / wireSegmentUnit_);
                    Point<double> legalLoc = legalizeOneBuffer(rawLoc, ffBufMaster);
                    bx = static_cast<int>(legalLoc.getX() * wireSegmentUnit_);
                    by = static_cast<int>(legalLoc.getY() * wireSegmentUnit_);
                  }
                  const std::string bName = "clkbuf_grp_"
                      + std::to_string(myGrpId) + "_s"
                      + std::to_string(numSinks) + "_"
                      + std::to_string(d);
                  ClockInst& grpBuf = clock_.addClockBuffer(
                      bName, ffBufMaster, bx, by);
                  cts3dDb_->setClockInstTier(grpBuf, ffTier);
                  addTreeLevelBuffer(&grpBuf);
                  cur->addInst(grpBuf);
                  ClockSubNet& next = clock_.addSubNet(
                      "clknet_grp_" + std::to_string(myGrpId) + "_s"
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
          logger_->info(CTS, 493,
              "V50 grouped branch {:2d}: grouped={} FFs ({} bufs), "
              "direct={}, holdSkip={}, branchMean={:.4f}ns uniform={}",
              idx, branchGroupedFFs, branchGroupBufs,
              branchDirectFFs, branchHoldSkip, absTgt,
              enableClusterUniform_ ? 1 : 0);
          return;  // skip legacy per-branch N-buffer code
        }

        // JYJ (2026-02-27) V41-fix: Per-FF relay buffer path (rewrote V39).
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

          // V41: Helper lambda to process a single FF with relay decision
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
            const double ffDelta = ffTarget - globalMeanTarget_;
            // V47: Pre-compute dx/dy for Elmore model (moved before nRelay decision)
            const int dx = ffInst->getX() - driverInst->getX();
            const int dy = ffInst->getY() - driverInst->getY();

            // V47: Elmore x_useful relay positioning.
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

            // Hold safety: check per-FF hold budget
            if (nRelay > 0) {
              std::string ffName = ffInst->getName();
              auto slashPos = ffName.rfind('/');
              if (slashPos != std::string::npos) {
                ffName = ffName.substr(0, slashPos);
              }
              auto hIt = perFfHoldBudget_.find(ffName);
              if (hIt != perFfHoldBudget_.end()) {
                while (nRelay > 0
                    && hIt->second - nRelay * perFfBufDelay_ < perFfHoldGuard_) {
                  --nRelay;
                }
                if (nRelay == 0) ++branchHoldSkip;
              }
            }

            // Debug stats (use delta, not absolute, for accurate coverage)
            totalLpDelay_debug += std::max(0.0, ffDelta);
            ++totalBranchSinks_debug;

            if (nRelay > 0) {
              // V47: dx/dy pre-computed above (before nRelay decision).
              // Single relay uses Elmore relayFrac; multi-relay uses equal spacing.
              // JYJ (2026-03-01) V47a: Get FF's tier for relay tier assignment.
              // Relay must be on the SAME tier as its target FF to avoid cross-tier
              // relay→FF nets. Cross-tier relay→FF nets cause GRT-0183 heap underflow
              // in 3D maze routing (observed in V43/V46/V47 for minority-tier FFs in
              // upper-dominant clusters). Using branch_tier (cluster dominant tier)
              // was wrong when FF is on the minority tier of the cluster.
              const int ff_tier_relay = cts3dDb_->getInstTier(ffInst->getDbInst());
              const int relay_tier = (ff_tier_relay >= 0) ? ff_tier_relay : branch_tier;

              ClockSubNet* cur = driverSubNet;
              for (int r = 0; r < nRelay; ++r) {
                // V47: 1-relay → Elmore x_useful fraction; else equal spacing
                const double frac = (nRelay == 1)
                    ? relayFrac
                    : static_cast<double>(r + 1) / (nRelay + 1);
                const int rx = driverInst->getX()
                    + static_cast<int>(dx * frac);
                const int ry = driverInst->getY()
                    + static_cast<int>(dy * frac);

                const std::string rName = "clkbuf_relay_"
                    + std::to_string(idx) + "_"
                    + std::to_string(numSinks) + "_"
                    + std::to_string(r);
                ClockInst& relay = clock_.addClockBuffer(
                    rName, bufMaster, rx, ry);
                cts3dDb_->setClockInstTier(relay, relay_tier);
                addTreeLevelBuffer(&relay);
                cur->addInst(relay);
                ClockSubNet& next = clock_.addSubNet(
                    "clknet_relay_" + std::to_string(idx) + "_"
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
            auto sinkIt = mapLocationToSink_.find(loc);
            if (sinkIt == mapLocationToSink_.end()) {
              constexpr double eps = 1e-6;
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
                  "V39 relay: Sink not found at ({:.10f}, {:.10f})",
                  loc.getX(), loc.getY());
            }
            ClockInst* sinkInst = sinkIt->second;

            // V41-fix BUG#1: Check if this is a cluster buffer
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
              "V39 relay branch {:2d}: relay={} FFs ({} bufs), "
              "direct={}, holdSkip={}, branchMean={:.4f}ns",
              idx, branchRelayFFs, branchRelayBufs,
              branchDirectFFs, branchHoldSkip, absTgt);
          return;  // skip legacy per-branch N-buffer code
        }

        // V32c-fix: use delta relative to globalMean to cap inter-branch differential
        const double delta = absTgt - globalMeanTarget_;
        int num_bufs = 0;
        if (delta > 0.0 && singleBufDelay_ > 1e-9) {
          num_bufs = static_cast<int>(std::round(delta / singleBufDelay_));
          num_bufs = std::max(0, std::min(num_bufs, maxLeafDelayBufs_));
        }

        // V32c-debug: accumulate LP-CTS gap stats (CTS-408 summary after loop)
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
              const std::string bufName = "clkbuf_delay_"
                                          + std::to_string(idx) + "_"
                                          + std::to_string(b);
              ClockInst& buf = clock_.addClockBuffer(bufName, delayBufMaster, bx, by);
              cts3dDb_->setClockInstTier(buf, branch_tier);
              addTreeLevelBuffer(&buf);

              // buf is SINK of curSub (curSub already has a driver)
              curSub->addInst(buf);

              // Create new subnet; first addInst sets buf as its driver
              ClockSubNet& nextSub = clock_.addSubNet(
                  "clknet_delay_" + std::to_string(idx) + "_" + std::to_string(b));
              nextSub.addInst(buf);  // buf becomes driver of nextSub

              curSub = &nextSub;
            }

            // Last subnet is the leaf subnet; connect FFs as sinks
            curSub->setLeafLevel(true);
            const std::vector<Point<double>>& sinkLocs
                = leafTopology.getBranchSinksLocations(idx);
            for (const Point<double>& loc : sinkLocs) {
              auto sinkIt = mapLocationToSink_.find(loc);
              // V37-fix: epsilon fallback for float precision mismatch
              if (sinkIt == mapLocationToSink_.end()) {
                constexpr double eps = 1e-6;
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
                logger_->error(CTS, 404, "Sink not found (delay buffer branch).");
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

        // V32c-debug CTS-406: branch has LP target but delta ≤ 0 → skipped by delta trigger.
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
          // V37-fix: float/double precision mismatch from multi-level
          // computeBranchSinks() double→float→double roundtrip.
          // If exact lookup fails, search for nearest entry within epsilon.
          if (sinkIt == mapLocationToSink_.end()) {
            constexpr double eps = 1e-6;
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
            logger_->error(CTS, 80,
                           "Sink not found at ({:.10f}, {:.10f}), mapSize={}",
                           loc.getX(), loc.getY(), mapLocationToSink_.size());
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

  // V32c-debug CTS-408: LP-CTS translation efficiency summary.
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

  // JYJ (2026-02-26) V39: Per-FF relay summary log
  if (enablePerFfRelay_ && totalBranchSinks_debug > 0) {
    const double cov = totalLpDelay_debug > 1e-12
        ? 100.0 * totalDelivered_debug / totalLpDelay_debug : 0.0;
    logger_->info(CTS, 426,
        "V39 per-FF relay summary: {} relay bufs across {}/{} branches, "
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
      "clkbuf_0", options_->getRootBuffer(), centerX, centerY);

  // clang-format off
  if (center != legalCenter) {
    debugPrint(logger_, CTS, "legalizer", 2, "createSingleBufferClockNet "
	       "legalizeOneBuffer clkbuf_0: {} => {}", center, legalCenter);
  }
  // clang-format on

  addTreeLevelBuffer(&rootBuffer);
  ClockSubNet& clockSubNet = clock_.addSubNet("clknet_0");
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
      // JYJ (2026-02-06) Replaced mapBufferMasterToTier with Cts3DDatabase
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
      // JYJ (2026-02-07) Tier assignment through SSOT (Cts3DDatabase)
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

  // JYJ (2026-02-06) Replaced mapBufferMasterToTier with Cts3DDatabase
  const std::string tieredMaster
      = tree_->getCts3DDatabase()->getBufferForTier(master, targetTier_);
  ClockInst& newBuffer
      = clock_->addClockBuffer(instPrefix_ + "_f",
                               tieredMaster,
                               target_.getX() * techCharDistUnit_,
                               target_.getY() * techCharDistUnit_);
  // JYJ (2026-02-07) Tier assignment through SSOT (Cts3DDatabase)
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

}  // namespace cts
