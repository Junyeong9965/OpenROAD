// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#pragma once

#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "CtsObserver.h"
#include "CtsOptions.h"
#include "TreeBuilder.h"
#include "Util.h"
#include "odb/db.h"
#include "odb/isotropy.h"

namespace cts {
class Graphics;

class SegmentBuilder
{
 public:
  SegmentBuilder(const std::string& instPrefix,
                 const std::string& netPrefix,
                 const Point<double>& root,
                 const Point<double>& target,
                 const std::vector<unsigned>& techCharWires,
                 Clock& clock,
                 ClockSubNet& drivingSubNet,
                 const TechChar& techChar,
                 unsigned techCharDistUnit,
                 TreeBuilder* tree,
                 odb::dbDatabase* db,
                 int targetTier);

  void build(const std::string& forceBuffer = "");
  void forceBufferInSegment(const std::string& master);

  ClockSubNet* getDrivingSubNet() const { return drivingSubNet_; }
  unsigned getNumBufferLevels() const { return numBufferLevels_; }
  TreeBuilder* getTree() const { return tree_; }

 private:
  const std::string instPrefix_;
  const std::string netPrefix_;
  const Point<double> root_;
  const Point<double> target_;
  const std::vector<unsigned> techCharWires_;
  const TechChar* techChar_;
  const unsigned techCharDistUnit_;
  Clock* clock_;
  ClockSubNet* drivingSubNet_;
  TreeBuilder* tree_;
  odb::dbDatabase* db_;
  int targetTier_ = -1;  // JYJ (2026-02-06) Tier for buffer selection in 3D CTS
  unsigned numBufferLevels_ = 0;
};

//-----------------------------------------------------------------------------
class HTreeBuilder : public TreeBuilder
{
  class LevelTopology
  {
   public:
    static constexpr unsigned NO_PARENT = std::numeric_limits<unsigned>::max();

    explicit LevelTopology(double length) : length_(length) {}

    void addWireSegment(unsigned idx) { wireSegments_.push_back(idx); }

    unsigned addBranchingPoint(const Point<double>& loc, unsigned parent)
    {
      branchPointLoc_.push_back(loc);
      parents_.push_back(parent);
      branchSinkLocs_.resize(branchPointLoc_.size());
      branchDrivingSubNet_.resize(branchPointLoc_.size(), nullptr);
      return branchPointLoc_.size() - 1;
    }

    void addSinkToBranch(unsigned branchIdx, const Point<double>& sinkLoc)
    {
      branchSinkLocs_[branchIdx].push_back(sinkLoc);
    }

    unsigned getBranchingPointSize() { return branchPointLoc_.size(); }

    Point<double>& getBranchingPoint(unsigned idx)
    {
      return branchPointLoc_[idx];
    }

    unsigned getBranchingPointParentIdx(unsigned idx) const
    {
      return parents_[idx];
    }

    double getLength() const { return length_; }
    void setLength(double x) { length_ = x; }

    void forEachBranchingPoint(
        const std::function<void(unsigned, Point<double>)>& func) const
    {
      for (unsigned idx = 0; idx < branchPointLoc_.size(); ++idx) {
        func(idx, branchPointLoc_[idx]);
      }
    }

    ClockSubNet* getBranchDrivingSubNet(unsigned idx) const
    {
      return branchDrivingSubNet_[idx];
    }

    void setBranchDrivingSubNet(unsigned idx, ClockSubNet& subNet)
    {
      branchDrivingSubNet_[idx] = &subNet;
    }

    const std::vector<unsigned>& getWireSegments() const
    {
      return wireSegments_;
    }

    const std::vector<Point<double>>& getBranchSinksLocations(
        unsigned branchIdx) const
    {
      return branchSinkLocs_[branchIdx];
    }

    void setOutputSlew(unsigned slew) { outputSlew_ = slew; }
    unsigned getOutputSlew() const { return outputSlew_; }
    void setOutputCap(unsigned cap) { outputCap_ = cap; }
    unsigned getOutputCap() const { return outputCap_; }
    void setRemainingLength(unsigned length) { remainingLength_ = length; }
    unsigned getRemainingLength() const { return remainingLength_; }
    void setCurrWl(int wl) { curr_Wl_ = wl; }
    int getCurrWl() const { return curr_Wl_; }

   private:
    double length_;
    unsigned outputSlew_ = 0;
    unsigned outputCap_ = 0;
    unsigned remainingLength_ = 0;
    int curr_Wl_ = 0;
    std::vector<unsigned> wireSegments_;
    std::vector<Point<double>> branchPointLoc_;
    std::vector<unsigned> parents_;
    std::vector<ClockSubNet*> branchDrivingSubNet_;
    std::vector<std::vector<Point<double>>> branchSinkLocs_;
  };

 public:
  HTreeBuilder(CtsOptions* options,
               Clock& net,
               TreeBuilder* parent,
               utl::Logger* logger,
               odb::dbDatabase* db)
      : TreeBuilder(options, net, parent, logger, db)
  {
  }

  void run() override;
  Point<double> legalizeOneBuffer(Point<double> bufferLoc,
                                  const std::string& bufferName) override;
  void findLegalLocations(const Point<double>& parentPoint,
                          const Point<double>& branchPoint,
                          double x1,
                          double y1,
                          double x2,
                          double y2,
                          std::vector<Point<double>>& points);
  void addCandidateLoc(double x,
                       double y,
                       const Point<double>& parentPoint,
                       double x1,
                       double y1,
                       double x2,
                       double y2,
                       std::vector<Point<double>>& points)
  {
    Point<double> candidate(x, y);
    if ((candidate != parentPoint) && isAlongBbox(x, y, x1, y1, x2, y2)) {
      points.emplace_back(x, y);
    }
  }
  void addCandidatePointsAlongBlockage(const Point<double>& point,
                                       const Point<double>& parentPoint,
                                       double targetDist,
                                       int scalingFactor,
                                       std::vector<Point<double>>& candidates,
                                       odb::Direction2D direction);
  Point<double> findBestLegalLocation(
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
      odb::Direction2D direction);
  Point<double> adjustBestLegalLocation(double targetDist,
                                        const Point<double>& currLoc,
                                        const Point<double>& parentPoint,
                                        const std::vector<Point<double>>& sinks,
                                        double x1,
                                        double y1,
                                        double x2,
                                        double y2,
                                        int scalingFactor,
                                        odb::Direction2D direction);
  void checkLegalityAndCostSpecial(const Point<double>& oldLoc,
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
                                   double& bestSinkDist);
  bool adjustAlongBlockage(double targetDist,
                           const Point<double>& currLoc,
                           const Point<double>& parentPoint,
                           const std::vector<Point<double>>& sinks,
                           double x1,
                           double y1,
                           double x2,
                           double y2,
                           int scalingFactor,
                           Point<double>& bestLoc);
  Point<double> adjustBeyondBlockage(const Point<double>& branchPoint,
                                     const Point<double>& parentPoint,
                                     double targetDist,
                                     const std::vector<Point<double>>& sinks,
                                     int scalingFactor,
                                     odb::Direction2D direction);
  void checkLegalityAndCost(const Point<double>& oldLoc,
                            const Point<double>& newLoc,
                            const Point<double>& parentPoint,
                            double targetDist,
                            const std::vector<Point<double>>& sinks,
                            int scalingFactor,
                            Point<double>& bestLoc,
                            double& sinkDist,
                            double& bestSinkDist);
  void legalize();
  void legalizeDummy();
  void printHTree();
  void plotSolution();
  std::string plotHTree();
  unsigned findSibling(LevelTopology& topology, unsigned i, unsigned par);
  Point<double>& findSiblingLoc(LevelTopology& topology,
                                unsigned i,
                                unsigned par)
  {
    unsigned j = findSibling(topology, i, par);
    return topology.getBranchingPoint(j);
  }

  std::vector<LevelTopology> getTopologyVector() const
  {
    return topologyForEachLevel_;
  }

  Box<double> getSinkRegion() const { return sinkRegion_; }

  int getWireSegmentUnit() const { return wireSegmentUnit_; }

  unsigned computeMinDelaySegment(unsigned length,
                                  unsigned inputSlew,
                                  unsigned inputCap,
                                  unsigned slewThreshold,
                                  int wirelengthThreshold,
                                  unsigned tolerance,
                                  unsigned& outputSlew,
                                  unsigned& outputCap,
                                  int& currWl) const;

  // JYJ (2026-02-09) Phase 3: Override computeDist to add HB penalty for cross-tier pairs
  double computeDist(const Point<double>& x, const Point<double>& y) override;

 private:
  void initSinkRegion();
  void computeLevelTopology(unsigned level, double width, double height);
  unsigned computeNumberOfSinksPerSubRegion(unsigned level) const;
  void computeSubRegionSize(unsigned level,
                            double& width,
                            double& height) const;
  unsigned computeMinDelaySegment(unsigned length) const;
  unsigned computeMinDelaySegment(unsigned length,
                                  unsigned inputSlew,
                                  unsigned inputCap,
                                  unsigned slewThreshold,
                                  unsigned tolerance,
                                  unsigned& outputSlew,
                                  unsigned& outputCap,
                                  bool forceBuffer,
                                  int expectedLength) const;
  void reportWireSegment(unsigned key) const;
  void createClockSubNets();
  void createSingleBufferClockNet();
  void initTopLevelSinks(std::vector<std::pair<float, float>>& sinkLocations,
                         std::vector<const ClockInst*>& sinkInsts);
  void initSecondLevelSinks(std::vector<std::pair<float, float>>& sinkLocations,
                            std::vector<const ClockInst*>& sinkInsts);
  void computeBranchSinks(
      const LevelTopology& topology,
      unsigned branchIdx,
      std::vector<std::pair<float, float>>& sinkLocations) const;

  bool isVertical(unsigned level) const
  {
    return level % 2 == (sinkRegion_.getHeight() >= sinkRegion_.getWidth());
  }
  bool isHorizontal(unsigned level) const { return !isVertical(level); }

  unsigned computeGridSizeX(unsigned level) const
  {
    return std::pow(2, (level + 1) / 2);
  }
  unsigned computeGridSizeY(unsigned level) const
  {
    return std::pow(2, level / 2);
  }

  void computeBranchingPoints(unsigned level, LevelTopology& topology);
  void refineBranchingPointsWithClustering(
      LevelTopology& topology,
      unsigned level,
      unsigned branchPtIdx1,
      unsigned branchPtIdx2,
      const Point<double>& rootLocation,
      const std::vector<std::pair<float, float>>& sinks);
  void preClusteringOpt(const std::vector<std::pair<float, float>>& sinks,
                        std::vector<std::pair<float, float>>& points,
                        std::vector<unsigned>& mapSinkToPoint);
  void preSinkClustering(const std::vector<std::pair<float, float>>& sinks,
                         const std::vector<const ClockInst*>& sinkInsts,
                         float maxDiameter,
                         unsigned clusterSize,
                         bool secondLevel = false);
  void assignSinksToBranches(
      LevelTopology& topology,
      unsigned branchPtIdx1,
      unsigned branchPtIdx2,
      const std::vector<std::pair<float, float>>& sinks,
      const std::vector<std::pair<float, float>>& points,
      const std::vector<unsigned>& mapSinkToPoint,
      const std::vector<std::vector<unsigned>>& clusters);

  bool isSubRegionTooSmall(double width, double height) const
  {
    return width < minLengthSinkRegion_ || height < minLengthSinkRegion_;
  }

  bool isNumberOfSinksTooSmall(unsigned numSinksPerSubRegion) const;

  double weightedDistance(const Point<double>& newLoc,
                          const Point<double>& oldLoc,
                          const std::vector<Point<double>>& sinks);
  void scalePosition(Point<double>& loc,
                     const Point<double>& parLoc,
                     double leng,
                     double scale);
  void adjustToplevelTopology(Point<double>& a,
                              Point<double>& b,
                              const Point<double>& parLoc);
  std::vector<unsigned> clusterDiameters() const { return clusterDiameters_; }
  std::vector<unsigned> clusterSizes() const { return clusterSizes_; }
  Point<double> resolveLocationCollision(
      const Point<double>& legalCenter) const;
  // JYJ (2026-02-06) Removed getDominantTierFromInsts, getDominantTierFromSinkLocs,
  // getDominantTierFromClockSinks - moved to Cts3DDatabase

 private:
  Box<double> sinkRegion_;
  std::vector<LevelTopology> topologyForEachLevel_;
  std::map<Point<double>, ClockInst*> mapLocationToSink_;
  std::vector<std::pair<float, float>> topLevelSinksClustered_;

  int wireSegmentUnit_ = 0;
  unsigned minInputCap_ = 0;
  unsigned numMaxLeafSinks_ = 0;
  unsigned minLengthSinkRegion_ = 0;
  unsigned clockTreeMaxDepth_ = 0;
  static constexpr int min_clustering_sinks_ = 200;
  static constexpr int min_clustering_macro_sinks_ = 10;
  std::vector<unsigned> clusterDiameters_ = {50, 100, 200};
  std::vector<unsigned> clusterSizes_ = {10, 20, 30};

  // JYJ (2026-02-21) Step 3: Skew-aware clustering penalty weight
  // β · |target_a - target_b| is added to computeDist() to group FFs
  // with similar arrival targets into the same cluster.
  double skewTargetBeta_ = 0.0;

  // JYJ (2026-02-21) Step 4: Per-branch delay targets
  // Maps (levelIdx, branchIdx) → mean arrival target for that branch.
  // Used in createClockSubNets() to insert/skip delay buffers.
  std::map<std::pair<int, unsigned>, double> branchDelayTargets_;
  double globalMeanTarget_ = 0.0;
  double delayTargetThreshold_ = 0.0;  // ns, from CTS_DELAY_TARGET_THRESHOLD env
  void computeBranchDelayTargets();

  // JYJ (2026-02-23) V31: Useful-skew wire length adjustment
  // Instead of a balanced (equal-wire) H-tree, shift each branch point
  // radially from the clock root based on the cluster's mean LP skew target.
  //   shift_um = wireSkewScale_ * (T_cluster - T_global) / wireDelayPerUnit_
  // Positive T_cluster → branch moves farther from root → longer wire → later clock.
  // Negative T_cluster → branch moves closer to root → shorter wire → earlier clock.
  // CTS_WIRE_SKEW_SCALE env var (default 1.0, 0 = off / zero-skew mode)
  // CTS_WIRE_DELAY_PS_UM env var: wire+buffer delay per unit length (ps/um, default 1.0)
  double wireSkewScale_    = 0.0;  // disabled by default (backward compatible)
  double wireDelayPerUnit_ = 0.001; // ns/um = 1.0 ps/um (ASAP7 typical CTS layer)

  // JYJ (2026-02-23) V32b: N-buffer chain leaf delay parameters.
  // Step 4 inserts N = round(branchMeanTarget / singleBufDelay_) delay buffers
  // in series at each leaf branch, instead of the V32 binary (0 or 1) approach.
  // Driven by absolute LP-TNS target (not relative to globalMean).
  //   CTS_LEAF_BUF_DELAY_NS  : delay per single buffer cell in ns (default 12ps)
  //   CTS_MAX_LEAF_DELAY_BUFS: max buffers per leaf branch (default 3)
  double singleBufDelay_   = 0.012;  // ns, from CTS_LEAF_BUF_DELAY_NS
  int    maxLeafDelayBufs_ = 3;      // from CTS_MAX_LEAF_DELAY_BUFS

  // JYJ (2026-02-23) V32c: cluster buffer → mean LP skew target map.
  // preSinkClustering() replaces individual FF positions in mapLocationToSink_
  // with cluster-buffer ClockInsts (clkbuf_leaf_X). These buffer names are NOT
  // in the LP CSV, so getSkewTarget(bufName) always returns 0.
  // Fix: compute mean FF target per cluster at clustering time and store here.
  // getClockInstTarget() uses this map for cluster buffers, LP CSV for raw FFs.
  std::unordered_map<ClockInst*, double> clusterBufTarget_;
  double getClockInstTarget(ClockInst* inst) const;

  // JYJ (2026-02-25) V37: Target-aware adaptive sub-clustering parameters.
  // After SinkClustering groups FFs spatially, clusters with high LP target
  // spread are split into sub-clusters so that per-leaf N-chain can cover
  // the reduced spread effectively.
  //   spread > split2way → 2-way split at median
  //   spread > split3way → 3-way split at terciles
  // CTS_ENABLE_TARGET_SPLIT (default 1): enable/disable
  // CTS_SPLIT_THRESHOLD_2WAY_PS (default 20): 2-way threshold in ps
  // CTS_SPLIT_THRESHOLD_3WAY_PS (default 40): 3-way threshold in ps
  bool   enableTargetSplit_      = true;
  double splitThreshold2wayNs_   = 0.020;  // ns
  double splitThreshold3wayNs_   = 0.040;  // ns

  // JYJ (2026-03-01) V48: Tier-aware cluster split.
  // If a leaf cluster has both bottom and upper FFs, split into two sub-clusters
  // so each gets the correct tier buffer (BUF_X4_bottom / BUF_X4_upper).
  // CTS_ENABLE_TIER_SPLIT (default 0): enable/disable
  bool enableTierSplit_ = false;

  // JYJ (2026-02-25) V37: Per-level globalMean for intermediate delay buffers.
  // computeBranchDelayTargets() now computes targets at ALL levels (not just leaf).
  // perLevelGlobalMean_[levelIdx] = weighted average of branch means at that level.
  std::map<int, double> perLevelGlobalMean_;

  // JYJ (2026-02-26) V39: Per-FF relay buffer parameters.
  // Replace per-BRANCH N-buffer (coverage=10.9%) with per-FF relay chains.
  // Leaf buffer -> [relay_0] -> [relay_1] -> FF  (high-target FFs)
  // Leaf buffer -> FF                            (low/zero-target FFs)
  // H-tree structure (62 leaf buffers) preserved; relays added after leaf only.
  //   CTS_ENABLE_PER_FF_RELAY:    0=off (legacy per-branch), 1=on
  //   CTS_PER_FF_BUF_DELAY_NS:   actual Liberty delay per relay buffer (ns)
  //   CTS_PER_FF_MAX_RELAY:       max relay buffers per FF
  //   CTS_PER_FF_HOLD_GUARD_NS:  min residual hold slack after relay (ns)
  //   CTS_PER_FF_MIN_TARGET_NS:   skip FFs with delta < this (ns)
  //   CTS_TIMING_GRAPH_CSV:       path to ff_timing_graph.csv for hold budgets
  // JYJ (2026-02-28) V47: Elmore-based x_useful relay positioning.
  // Replaces equal interpolation with exact position x_useful:
  //   delta(x) = d_buf - rc * x * (L - x)   [ps, rc = rw[kOhm/um]*cw[fF/um]]
  //   x_useful  = [L - sqrt(L^2 - 4*(d_buf-target)/rc)] / 2  [um]
  // Single relay at x_useful delivers target_delta exactly (no quantization).
  // Falls back to multi-relay equal spacing when target > d_buf.
  //   CTS_RELAY_RW_PER_UM:  wire resistance (kOhm/um); 0 = Elmore disabled
  //   CTS_RELAY_CW_PER_UM:  wire capacitance (fF/um)
  //   CTS_DBU_PER_UM:       database units per micron (default 2000)
  bool   enablePerFfRelay_     = false;
  double perFfBufDelay_        = 0.015;   // ns (V47: corrected to actual 15ps)
  int    perFfMaxRelay_        = 3;
  double perFfHoldGuard_       = 0.020;   // ns
  double perFfMinTarget_       = 0.010;   // ns
  // V47: Elmore wire RC parameters for x_useful computation
  double relayRwKOhmPerUm_     = 0.0;    // kOhm/um (0 = Elmore disabled)
  double relayCwFfPerUm_       = 0.0;    // fF/um
  double relayDbuPerUm_        = 2000.0; // DB units per micron

  // Per-FF hold budget: ff_name -> min hold slack (ns) as capture FF
  std::unordered_map<std::string, double> perFfHoldBudget_;
  void loadPerFfHoldBudgets(const std::string& timingGraphCsv);
  // JYJ (2026-02-26) V40: Load IO edge hold budgets (PI→FF hold slack)
  // to prevent relay insertion on FFs with tight PI→FF hold.
  void loadIoHoldBudgets(const std::string& ioCsv);

  // JYJ (2026-03-05) V49: Grouped Delay Chain parameters.
  // Instead of per-FF relay (1 chain per FF → many nets → GRT-0183),
  // group FFs within each leaf cluster by quantized LP target delta.
  // Group 0 (delta < min_delta): direct connection to leaf buffer.
  // Group k (k > 0): shared chain of k delay buffers → all group FFs.
  //   CTS_ENABLE_GROUPED_DELAY:         0/1 (default 0)
  //   CTS_GROUPED_DELAY_MAX_DEPTH:      max chain depth (default 3)
  //   CTS_GROUPED_DELAY_MIN_DELTA_NS:   min delta to create group (default 0.010)
  bool   enableGroupedDelay_       = false;
  int    groupedDelayMaxDepth_     = 3;
  double groupedDelayMinDelta_     = 0.010;  // ns (10ps)

  // JYJ (2026-03-06) V50: Cluster-Uniform Depth mode.
  // When enabled, all FFs in a leaf cluster share a SINGLE chain depth k_c
  // computed from the cluster median LP target delta.
  // V49a created up to MAX_DEPTH distinct depth groups per cluster →
  // 155 extra buffers (62 clusters × ~2.5 avg groups).
  // V50 creates exactly 1 chain per cluster (depth=k_c or 0) →
  // ~62 extra buffers max (1 chain per cluster), reducing buffer count ~60%.
  // Hold safety: k_c capped by min(hold_budget) across all FFs in cluster.
  //   CTS_GROUPED_DELAY_CLUSTER_UNIFORM: 0/1 (default 0; set 1 for V50 mode)
  bool enableClusterUniform_       = false;

  // V39: Env guard for V37 intermediate delay buffers (default off)
  bool enableMidDelayBufs_ = false;  // CTS_ENABLE_MID_DELAY_BUFS

  // JYJ (2026-02-26) V40: CKMeans target-aware branching weight.
  // When > 0, CKMeans distance includes target penalty so H-tree branching
  // groups FFs with similar LP targets together → better V31 wire shift.
  double ckmeansTargetBeta_ = 0.0;  // CTS_CKMEANS_TARGET_BETA
};

}  // namespace cts
