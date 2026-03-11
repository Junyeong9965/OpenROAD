// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include <cmath>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "utl/Logger.h"

namespace cts::CKMeans {

struct Sink;

class Clustering
{
 public:
  Clustering(const std::vector<std::pair<float, float>>& sinks,
             utl::Logger* logger);
  Clustering(const std::vector<std::pair<float, float>>& sinks,
             float xBranch,
             float yBranch,
             utl::Logger* logger);
  ~Clustering();

  void iterKmeans(unsigned iter,
                  unsigned n,
                  unsigned cap,
                  unsigned max,
                  unsigned power,
                  std::vector<std::pair<float, float>>& means);

  void getClusters(std::vector<std::vector<unsigned>>& newClusters) const;

  // JYJ (2026-02-26) V40: Set per-sink LP targets and distance weight.
  // When beta > 0, calcDist adds target penalty to cluster assignments.
  void setSinkTargets(const std::vector<float>& targets, float beta);

 private:
  float Kmeans(unsigned n,
               unsigned cap,
               unsigned max,
               unsigned power,
               std::vector<std::pair<float, float>>& means);
  float calcSilh(const std::vector<std::pair<float, float>>& means) const;
  void minCostFlow(const std::vector<std::pair<float, float>>& means,
                   unsigned cap,
                   float dist,
                   unsigned power);
  void fixSegmentLengths(std::vector<std::pair<float, float>>& means);
  void fixSegment(const std::pair<float, float>& fixedPoint,
                  float targetDist,
                  std::pair<float, float>& movablePoint);

  // JYJ V40: target-aware distance (non-static, uses targetBeta_ + meanTargets_)
  float calcDist(const std::pair<float, float>& loc, size_t clusterIdx,
                 const Sink* sink) const;
  // Pure geometric distance (static, for fixSegment wire normalization)
  static float calcDist(const std::pair<float, float>& loc1,
                        const std::pair<float, float>& loc2);

  utl::Logger* logger_;
  std::vector<Sink> sinks_;
  std::vector<std::vector<Sink*>> clusters_;

  float segment_length_ = 0.0;
  std::optional<std::pair<float, float>> branching_point_;

  // JYJ V40: LP target-aware clustering parameters
  float targetBeta_{0.0f};          // target penalty weight in distance
  std::vector<float> meanTargets_;  // per-cluster mean target (updated each iter)
};

}  // namespace cts::CKMeans
