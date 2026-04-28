// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
//
// TARP: FM min-cut with Hilbert initial partition
// + degree-0 pinning.
//
// Algorithm overview:
//   1. Build timing affinity graph (slack-filtered, top-K pruned)
//   2. Compute global Hilbert indices (one-time, O(N log N))
//   3. Recursive FM bisection:
//      a. Initial partition by Hilbert ordering (capacity-balanced)
//      b. Pin degree-0 nodes in current subgraph (per-recursion-level)
//      c. FM passes on unpinned nodes only
//      d. Base case: size <= max_cluster_size AND diameter <= maxDiameter
//   4. S_max validation on leaf clusters

#include "TarpClustering.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

namespace cts {

TarpClustering::TarpClustering(utl::Logger* logger, const TarpParams& params)
    : logger_(logger), params_(params)
{
  s_max_ = (params_.max_depth + 1) * params_.d_buf;
}

void TarpClustering::addSink(unsigned orig_idx, double x, double y, float cap,
                              double lp_target, double hold_budget,
                              const std::string& name)
{
  sinks_.push_back({orig_idx, x, y, cap, lp_target, hold_budget, name});
}

void TarpClustering::addTimingEdge(unsigned ff_a, unsigned ff_b,
                                    double setup_slack, double T_period)
{
  if (ff_a >= sinks_.size() || ff_b >= sinks_.size() || ff_a == ff_b) return;

  // Slack threshold filter: skip non-critical edges.
  // Dense designs (IBEX 153 edges/FF) -> most edges pruned -> FM focuses on critical pairs.
  // Sparse designs (AES 15 edges/FF) -> most edges kept -> full timing awareness.
  const double slack_thresh = params_.slack_threshold_ratio * T_period;
  if (setup_slack > slack_thresh) return;

  // Lazy init on first edge
  if (adj_.size() < sinks_.size()) {
    adj_.resize(sinks_.size());
    computeDieDiagonal();
  }

  const double dx = std::abs(sinks_[ff_a].x - sinks_[ff_b].x);
  const double dy = std::abs(sinks_[ff_a].y - sinks_[ff_b].y);
  const double w_spatial = std::max(0.0, 1.0 - (dx + dy) / die_diagonal_);

  // Continuous slack-proportional weight:
  // more negative slack → higher weight (up to 2.0)
  // removes binary jump at s=0 that caused
  // over-aggressive FM moves for small violations
  double w_timing = std::min(
      2.0, std::max(0.0, 1.0 - setup_slack / T_period));

  const double affinity = params_.alpha * w_spatial
                         + (1.0 - params_.alpha) * w_timing;

  adj_[ff_a].emplace_back(ff_b, affinity);
  adj_[ff_b].emplace_back(ff_a, affinity);
}

void TarpClustering::computeDieDiagonal()
{
  if (sinks_.empty()) { die_diagonal_ = 1.0; return; }
  double x_min = sinks_[0].x, x_max = x_min;
  double y_min = sinks_[0].y, y_max = y_min;
  for (const auto& s : sinks_) {
    x_min = std::min(x_min, s.x); x_max = std::max(x_max, s.x);
    y_min = std::min(y_min, s.y); y_max = std::max(y_max, s.y);
  }
  die_diagonal_ = std::sqrt((x_max-x_min)*(x_max-x_min)
                            + (y_max-y_min)*(y_max-y_min));
  if (die_diagonal_ < 1e-9) die_diagonal_ = 1.0;
}

void TarpClustering::buildAffinityGraph()
{
  // Top-K pruning per FF: keep highest affinity edges only.
  // No KNN spatial edges -- degree-0 FFs are pinned by Hilbert ordering.
  const int K = params_.top_k_edges;
  for (auto& neighbors : adj_) {
    if (static_cast<int>(neighbors.size()) > K) {
      std::partial_sort(neighbors.begin(), neighbors.begin() + K, neighbors.end(),
                        [](const auto& a, const auto& b) {
                          return a.second > b.second;
                        });
      neighbors.resize(K);
    }
  }

  size_t total = 0;
  unsigned degree0 = 0;
  for (const auto& a : adj_) {
    total += a.size();
    if (a.empty()) degree0++;
  }
  logger_->info(utl::CTS, 830,
      "TARP: {} FFs, {} edges (slack<{:.1f}*T, top-{}), {} degree-0 (pinned), "
      "diag={:.2f}",
      sinks_.size(), total, params_.slack_threshold_ratio, K,
      degree0, die_diagonal_);
}

// Standard Hilbert curve xy-to-index (Wikipedia implementation).
// n must be power of 2. x, y passed by value (mutated internally during rotation).
unsigned TarpClustering::xy2d(unsigned n, unsigned x, unsigned y)
{
  unsigned d = 0;
  for (unsigned s = n / 2; s > 0; s /= 2) {
    unsigned rx = (x & s) > 0 ? 1 : 0;
    unsigned ry = (y & s) > 0 ? 1 : 0;
    d += s * s * ((3 * rx) ^ ry);
    if (ry == 0) {
      if (rx == 1) { x = s - 1 - x; y = s - 1 - y; }
      std::swap(x, y);
    }
  }
  return d;
}

void TarpClustering::computeHilbertIndices()
{
  const unsigned N = sinks_.size();
  hilbert_idx_.resize(N);
  if (N == 0) return;

  double x_min = sinks_[0].x, x_max = x_min;
  double y_min = sinks_[0].y, y_max = y_min;
  for (const auto& s : sinks_) {
    x_min = std::min(x_min, s.x); x_max = std::max(x_max, s.x);
    y_min = std::min(y_min, s.y); y_max = std::max(y_max, s.y);
  }

  const unsigned grid_size = 1024;  // 2^10
  const double x_range = std::max(x_max - x_min, 1e-9);
  const double y_range = std::max(y_max - y_min, 1e-9);

  for (unsigned i = 0; i < N; ++i) {
    unsigned gx = static_cast<unsigned>(
        (sinks_[i].x - x_min) / x_range * (grid_size - 1));
    unsigned gy = static_cast<unsigned>(
        (sinks_[i].y - y_min) / y_range * (grid_size - 1));
    gx = std::min(gx, grid_size - 1);
    gy = std::min(gy, grid_size - 1);
    hilbert_idx_[i] = xy2d(grid_size, gx, gy);
  }
}

void TarpClustering::run()
{
  if (sinks_.empty()) return;
  adj_.resize(sinks_.size());
  computeDieDiagonal();
  buildAffinityGraph();
  computeHilbertIndices();

  std::vector<unsigned> all(sinks_.size());
  std::iota(all.begin(), all.end(), 0u);
  recursivePartition(all, 0);

  logger_->info(utl::CTS, 831, "TARP: {} clusters from {} FFs",
                clusters_.size(), sinks_.size());
}

void TarpClustering::recursivePartition(const std::vector<unsigned>& ffs,
                                         int depth)
{
  // Base case: small enough for a leaf cluster
  if (ffs.size() <= params_.max_cluster_size) {
    // maxDiameter guard: FM doesn't guarantee spatial compactness
    if (computeDiameter(ffs) > params_.max_diameter) {
      splitBySpatialMedian(ffs);
      return;
    }
    // S_max delivery range check
    if (validateCluster(ffs)) {
      clusters_.push_back(ffs);
    } else {
      splitByTarget(ffs);
    }
    return;
  }

  std::vector<unsigned> left, right;
  double cut = fmBisection(ffs, left, right);

  // Degenerate partition fallback: empty, or balance violation (cut == -1)
  if (left.empty() || right.empty() || cut < -0.5) {
    auto sorted = ffs;
    std::sort(sorted.begin(), sorted.end(), [this](unsigned a, unsigned b) {
      return hilbert_idx_[a] < hilbert_idx_[b];
    });
    size_t mid = sorted.size() / 2;
    left.assign(sorted.begin(), sorted.begin() + mid);
    right.assign(sorted.begin() + mid, sorted.end());
  }

  recursivePartition(left, depth + 1);
  recursivePartition(right, depth + 1);
}

double TarpClustering::fmBisection(const std::vector<unsigned>& ffs,
                                    std::vector<unsigned>& left,
                                    std::vector<unsigned>& right)
{
  const unsigned N = ffs.size();
  if (N <= 1) { left = ffs; return 0.0; }

  // Build local subgraph adjacency list
  auto local_adj = buildSubgraphAdjList(ffs);

  // --- Step 4a: Initial partition by Hilbert ordering ---
  std::vector<unsigned> sorted_local(N);
  std::iota(sorted_local.begin(), sorted_local.end(), 0u);
  std::sort(sorted_local.begin(), sorted_local.end(),
            [&](unsigned a, unsigned b) {
              return hilbert_idx_[ffs[a]] < hilbert_idx_[ffs[b]];
            });

  // Capacity-balanced split point (cumulative pin cap)
  float total_cap = 0;
  for (unsigned i = 0; i < N; ++i) total_cap += sinks_[ffs[i]].cap;
  const float target_cap = total_cap * 0.5f;
  const float cap_lo = total_cap * (1.0f - static_cast<float>(params_.balance_ratio));
  const float cap_hi = total_cap * static_cast<float>(params_.balance_ratio);

  std::vector<int> part(N, 1);  // default: partition B
  float cum_cap = 0;
  for (unsigned i = 0; i < N; ++i) {
    unsigned li = sorted_local[i];
    cum_cap += sinks_[ffs[li]].cap;
    if (cum_cap <= target_cap) {
      part[li] = 0;  // partition A
    }
  }

  // --- Step 4b: Degree-0 pinning (per-recursion-level subgraph) ---
  // FFs with no edges in THIS subgraph stay at Hilbert position.
  std::vector<bool> pinned(N, false);
  unsigned pin_count = 0;
  for (unsigned i = 0; i < N; ++i) {
    if (local_adj[i].empty()) {
      pinned[i] = true;
      pin_count++;
    }
  }

  // --- Step 4c: FM passes (unpinned nodes only) ---
  // Optimized: incremental cut tracking + incremental cap + O(degree) gain update.
  // Original was O(N*E) per pass due to computeCut() and sideCap() in inner loop.
  // Now O(E) per pass: cut delta computed from moved node's edges only.

  // Compute per-node gain: external - internal edge weight
  auto computeGain = [&](unsigned i) -> double {
    double internal = 0, external = 0;
    for (const auto& [nbr, wt] : local_adj[i]) {
      if (part[nbr] == part[i]) internal += wt; else external += wt;
    }
    return external - internal;
  };

  // Initial cut (one-time O(E) computation)
  double current_cut = 0;
  for (unsigned i = 0; i < N; ++i)
    for (const auto& [nbr, wt] : local_adj[i])
      if (part[nbr] != part[i]) current_cut += wt;
  current_cut /= 2.0;

  // Incremental cap tracking
  float cap_side0 = 0, cap_side1 = 0;
  for (unsigned i = 0; i < N; ++i) {
    if (part[i] == 0) cap_side0 += sinks_[ffs[i]].cap;
    else              cap_side1 += sinks_[ffs[i]].cap;
  }

  double best_cut = current_cut;
  auto best_part = part;

  for (int pass = 0; pass < params_.fm_max_passes; ++pass) {
    std::vector<double> gains(N);
    for (unsigned i = 0; i < N; ++i) {
      gains[i] = pinned[i] ? -1e30 : computeGain(i);
    }

    std::vector<bool> locked(N, false);
    for (unsigned i = 0; i < N; ++i) {
      if (pinned[i]) locked[i] = true;
    }

    // Reset incremental cap for this pass
    cap_side0 = 0; cap_side1 = 0;
    for (unsigned i = 0; i < N; ++i) {
      if (part[i] == 0) cap_side0 += sinks_[ffs[i]].cap;
      else              cap_side1 += sinks_[ffs[i]].cap;
    }

    // Recompute current_cut at pass start
    current_cut = 0;
    for (unsigned i = 0; i < N; ++i)
      for (const auto& [nbr, wt] : local_adj[i])
        if (part[nbr] != part[i]) current_cut += wt;
    current_cut /= 2.0;

    std::vector<std::pair<unsigned, int>> moves;
    double best_prefix_cut = current_cut;
    unsigned best_prefix = 0;

    for (unsigned step = 0; step < N - pin_count; ++step) {
      // Find best unlocked node with max gain (O(N) scan — acceptable for FM)
      int best_node = -1;
      double best_gain = -1e30;
      for (unsigned i = 0; i < N; ++i) {
        if (locked[i] || gains[i] <= best_gain) continue;
        // Incremental balance check: O(1) instead of O(N)
        float cap_i = sinks_[ffs[i]].cap;
        float from_cap = (part[i] == 0) ? cap_side0 : cap_side1;
        float to_cap   = (part[i] == 0) ? cap_side1 : cap_side0;
        if (from_cap - cap_i < cap_lo || to_cap + cap_i > cap_hi) continue;
        best_gain = gains[i];
        best_node = i;
      }
      if (best_node < 0) break;

      unsigned mv = static_cast<unsigned>(best_node);
      moves.emplace_back(mv, part[mv]);

      // Incremental cut update: moving mv from side S to side (1-S)
      // For each neighbor of mv: if same side → becomes cut edge (+wt)
      //                          if diff side → becomes internal (-wt)
      for (const auto& [nbr, wt] : local_adj[mv]) {
        if (part[nbr] == part[mv]) current_cut += wt;  // was internal, now cut
        else                       current_cut -= wt;  // was cut, now internal
      }

      // Incremental cap update
      float cap_mv = sinks_[ffs[mv]].cap;
      if (part[mv] == 0) { cap_side0 -= cap_mv; cap_side1 += cap_mv; }
      else                { cap_side1 -= cap_mv; cap_side0 += cap_mv; }

      part[mv] = 1 - part[mv];
      locked[mv] = true;

      // Update gains for unlocked neighbors only — O(degree) per neighbor
      for (const auto& [nbr, wt] : local_adj[mv]) {
        if (!locked[nbr]) gains[nbr] = computeGain(nbr);
      }

      if (current_cut < best_prefix_cut) {
        best_prefix_cut = current_cut;
        best_prefix = step + 1;
      }
    }

    // Rollback moves after best prefix + restore incremental state
    for (unsigned i = moves.size(); i > best_prefix; --i) {
      unsigned ri = moves[i-1].first;
      int old_side = moves[i-1].second;
      // Restore cut incrementally
      for (const auto& [nbr, wt] : local_adj[ri]) {
        if (part[nbr] == part[ri]) current_cut += wt;
        else                       current_cut -= wt;
      }
      // Restore cap
      float cap_ri = sinks_[ffs[ri]].cap;
      if (part[ri] == 0) { cap_side0 -= cap_ri; cap_side1 += cap_ri; }
      else                { cap_side1 -= cap_ri; cap_side0 += cap_ri; }
      part[ri] = old_side;
    }

    if (current_cut < best_cut - 1e-9) {
      best_cut = current_cut;
      best_part = part;
    } else {
      break;  // converged
    }
  }

  part = best_part;
  left.clear(); right.clear();
  for (unsigned i = 0; i < N; ++i) {
    if (part[i] == 0) left.push_back(ffs[i]); else right.push_back(ffs[i]);
  }

  // Balance violation check -> signal degenerate to caller
  float left_cap = 0, right_cap = 0;
  for (unsigned i : left) left_cap += sinks_[i].cap;
  for (unsigned i : right) right_cap += sinks_[i].cap;
  if (left_cap < cap_lo || right_cap < cap_lo
      || left_cap > cap_hi || right_cap > cap_hi) {
    left.clear(); right.clear();
    return -1.0;  // caller uses Hilbert median fallback
  }

  return best_cut;
}

bool TarpClustering::validateCluster(const std::vector<unsigned>& ffs) const
{
  if (ffs.size() <= 1) return true;
  double t_min = sinks_[ffs[0]].lp_target, t_max = t_min;
  for (unsigned idx : ffs) {
    t_min = std::min(t_min, sinks_[idx].lp_target);
    t_max = std::max(t_max, sinks_[idx].lp_target);
  }
  return (t_max - t_min) <= s_max_;
}

void TarpClustering::splitByTarget(const std::vector<unsigned>& ffs)
{
  if (ffs.size() <= 1) { clusters_.push_back(ffs); return; }
  auto sorted = ffs;
  std::sort(sorted.begin(), sorted.end(), [this](unsigned a, unsigned b) {
    return sinks_[a].lp_target < sinks_[b].lp_target;
  });
  size_t mid = sorted.size() / 2;
  std::vector<unsigned> lo(sorted.begin(), sorted.begin() + mid);
  std::vector<unsigned> hi(sorted.begin() + mid, sorted.end());

  auto processHalf = [this](const std::vector<unsigned>& half) {
    if (half.empty()) return;
    if (computeDiameter(half) > params_.max_diameter) {
      splitBySpatialMedian(half);
    } else if (validateCluster(half)) {
      clusters_.push_back(half);
    } else {
      splitByTarget(half);
    }
  };
  processHalf(lo);
  processHalf(hi);
}

void TarpClustering::splitBySpatialMedian(const std::vector<unsigned>& ffs)
{
  if (ffs.size() <= 1) { clusters_.push_back(ffs); return; }

  // Split along longer dimension
  double xn = sinks_[ffs[0]].x, xx = xn, yn = sinks_[ffs[0]].y, yx = yn;
  for (unsigned i : ffs) {
    xn = std::min(xn, sinks_[i].x); xx = std::max(xx, sinks_[i].x);
    yn = std::min(yn, sinks_[i].y); yx = std::max(yx, sinks_[i].y);
  }

  auto sorted = ffs;
  if ((xx - xn) >= (yx - yn)) {
    std::sort(sorted.begin(), sorted.end(),
              [this](unsigned a, unsigned b) { return sinks_[a].x < sinks_[b].x; });
  } else {
    std::sort(sorted.begin(), sorted.end(),
              [this](unsigned a, unsigned b) { return sinks_[a].y < sinks_[b].y; });
  }

  size_t mid = sorted.size() / 2;
  std::vector<unsigned> lo(sorted.begin(), sorted.begin() + mid);
  std::vector<unsigned> hi(sorted.begin() + mid, sorted.end());

  auto processHalf = [this](const std::vector<unsigned>& half) {
    if (half.empty()) return;
    if (half.size() <= params_.max_cluster_size
        && computeDiameter(half) <= params_.max_diameter
        && validateCluster(half)) {
      clusters_.push_back(half);
    } else if (half.size() <= params_.max_cluster_size) {
      if (computeDiameter(half) > params_.max_diameter) {
        splitBySpatialMedian(half);
      } else {
        splitByTarget(half);
      }
    } else {
      recursivePartition(half, 0);
    }
  };
  processHalf(lo);
  processHalf(hi);
}

double TarpClustering::computeDiameter(const std::vector<unsigned>& ffs) const
{
  if (ffs.size() <= 1) return 0.0;
  double xn = sinks_[ffs[0]].x, xx = xn, yn = sinks_[ffs[0]].y, yx = yn;
  for (unsigned i : ffs) {
    xn = std::min(xn, sinks_[i].x); xx = std::max(xx, sinks_[i].x);
    yn = std::min(yn, sinks_[i].y); yx = std::max(yx, sinks_[i].y);
  }
  return (xx - xn) + (yx - yn);
}

std::vector<std::vector<std::pair<unsigned, double>>>
TarpClustering::buildSubgraphAdjList(const std::vector<unsigned>& ffs) const
{
  const unsigned N = ffs.size();
  std::unordered_map<unsigned, unsigned> g2l;
  for (unsigned i = 0; i < N; ++i) g2l[ffs[i]] = i;

  std::vector<std::vector<std::pair<unsigned, double>>> la(N);
  for (unsigned i = 0; i < N; ++i) {
    for (const auto& [gnbr, wt] : adj_[ffs[i]]) {
      auto it = g2l.find(gnbr);
      if (it != g2l.end()) la[i].emplace_back(it->second, wt);
    }
  }
  return la;
}

}  // namespace cts
