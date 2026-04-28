// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
//
// TARP: FM min-cut with Hilbert initial partition
// + degree-0 pinning. Replaces KNN spatial edge approach.
//
// Key design decisions:
//   - Hilbert ordering as FM initial partition → spatial quality baseline
//   - Degree-0 nodes (no timing edges in subgraph) pinned → stay at Hilbert position
//   - No KNN spatial edges → timing signal not diluted
//   - Slack threshold filter → dense graphs auto-pruned
//   - maxDiameter guard on leaf clusters → spatial compactness guaranteed
//   - Per-recursion-level degree recomputation → correct pinning at each level

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "Util.h"
#include "utl/Logger.h"

namespace cts {

struct TarpEdge {
  unsigned ff_a;
  unsigned ff_b;
  double affinity;
};

struct TarpSink {
  unsigned orig_idx;
  double x, y;
  float cap;
  double lp_target;
  double hold_budget;
  std::string name;
};

struct TarpParams {
  double alpha = 0.7;
  unsigned max_cluster_size = 10;
  double max_diameter = 100.0;   // in wireSegmentUnit (not um)
  double d_buf = 0.012;
  int max_depth = 3;
  int top_k_edges = 30;
  double slack_threshold_ratio = 0.3;
  int fm_max_passes = 10;
  double balance_ratio = 0.6;
  int min_ffs_for_tarp = 50;
};

using TarpCluster = std::vector<unsigned>;

class TarpClustering
{
 public:
  TarpClustering(utl::Logger* logger, const TarpParams& params);

  void addSink(unsigned orig_idx, double x, double y, float cap,
               double lp_target, double hold_budget, const std::string& name);
  void addTimingEdge(unsigned ff_a, unsigned ff_b, double setup_slack,
                     double T_period);
  void run();

  const std::vector<TarpCluster>& clusters() const { return clusters_; }
  const std::vector<TarpSink>& sinks() const { return sinks_; }

 private:
  void buildAffinityGraph();
  void computeHilbertIndices();
  void recursivePartition(const std::vector<unsigned>& ff_indices, int depth);
  double fmBisection(const std::vector<unsigned>& ff_indices,
                     std::vector<unsigned>& left,
                     std::vector<unsigned>& right);
  bool validateCluster(const std::vector<unsigned>& cluster_ffs) const;
  void splitByTarget(const std::vector<unsigned>& cluster_ffs);
  void splitBySpatialMedian(const std::vector<unsigned>& ffs);
  double computeDiameter(const std::vector<unsigned>& ff_indices) const;
  void computeDieDiagonal();
  std::vector<std::vector<std::pair<unsigned, double>>>
  buildSubgraphAdjList(const std::vector<unsigned>& ff_indices) const;

  // Hilbert curve xy-to-index (standard Wikipedia implementation).
  // x, y are passed by value (mutated internally during rotation).
  static unsigned xy2d(unsigned n, unsigned x, unsigned y);

  utl::Logger* logger_;
  TarpParams params_;
  std::vector<TarpSink> sinks_;
  std::vector<std::vector<std::pair<unsigned, double>>> adj_;
  std::vector<TarpCluster> clusters_;
  std::vector<unsigned> hilbert_idx_;  // per-FF Hilbert curve index
  double die_diagonal_ = 1.0;
  double s_max_ = 0.0;
};

}  // namespace cts
