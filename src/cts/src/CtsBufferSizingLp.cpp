// SPDX-License-Identifier: BSD-3-Clause
// C++ buffer sizing LP using OR-Tools GLOP
// Ported from Python buffer_sizing_lp.py (V35 formulation)

#include "CtsBufferSizingLp.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <set>

#include "odb/db.h"
#include "ortools/linear_solver/linear_solver.h"
#include "sta/Liberty.hh"
#include "sta/Network.hh"
#include "sta/Sta.hh"
#include "utl/Logger.h"

namespace cts {

using operations_research::MPConstraint;
using operations_research::MPObjective;
using operations_research::MPSolver;
using operations_research::MPVariable;

// Log IDs for buffer sizing (CTS-610 series)
constexpr int CTS_BS_START = 610;
constexpr int CTS_BS_BUFS = 611;
constexpr int CTS_BS_LP = 612;
constexpr int CTS_BS_RESULT = 613;

CtsBufferSizingLp::CtsBufferSizingLp(odb::dbBlock* block,
                                     sta::dbSta* sta,
                                     sta::dbNetwork* network,
                                     utl::Logger* logger)
    : block_(block), sta_(sta), network_(network), logger_(logger)
{
}

bool CtsBufferSizingLp::isClockBuffer(odb::dbInst* inst) const
{
  // Clock buffers match patterns: clkbuf_*, clkload*, sg_*, *relay*
  std::string name = inst->getName();
  return (name.find("clkbuf_") == 0
       || name.find("clkload") == 0
       || name.find("sg_leaf_") == 0
       || name.find("sg_trunk_") == 0
       || name.find("sg_root") == 0
       || name.find("sg_dly_") == 0
       || name.find("relay") != std::string::npos
       || name.find("clkbuf_ghtap_") == 0);
}

float CtsBufferSizingLp::getBufferDelay(const std::string& master_name) const
{
  // Get buffer intrinsic delay from Liberty
  // Uses Liberty intrinsic delay
  // For now, use simple model: delay ∝ 1/sqrt(strength)
  // Calibrated: BUFx4 = 12ps, BUFx8 = 8.5ps, BUFx2 = 17ps
  std::string base = master_name.substr(0, master_name.find('_'));
  int strength = 4;
  try {
    std::string num_str = base.substr(4);  // "BUFx8" → "8"
    // Remove trailing 'f' if present
    if (!num_str.empty() && num_str.back() == 'f') {
      num_str.pop_back();
    }
    strength = std::stoi(num_str);
  } catch (...) {
    strength = 4;
  }
  // Model: delay = 24 / sqrt(strength) ps
  return 24.0f / std::sqrt(static_cast<float>(strength));
}

void CtsBufferSizingLp::extractTimingEdges()
{
  // Extract post-CTS FF→FF timing edges from STA
  // This uses the same approach as VerilogFFExtractor::fillTimingInfo
  // but operates on the post-CTS database directly
  logger_->info(utl::CTS, CTS_BS_START,
                "BufferSizing: Extracting post-CTS timing edges...");

  // Post-CTS edges loaded from CSV
  // For now, edges will be populated externally from Tcl
  logger_->info(utl::CTS, CTS_BS_START + 1,
                "BufferSizing: {} timing edges loaded", edges_.size());
}

void CtsBufferSizingLp::extractBufferInfo()
{
  // Find all clock buffers in the design
  buffers_.clear();
  ff_to_buf_.clear();

  for (auto* inst : block_->getInsts()) {
    if (!isClockBuffer(inst)) continue;

    ClockBufInfo buf;
    buf.inst_name = inst->getName();
    buf.cell_name = inst->getMaster()->getName();
    buf.tier = (buf.cell_name.find("_bottom") != std::string::npos) ? "bottom" : "upper";

    // Find driven FFs by tracing output net
    auto* oterm = inst->findITerm("Y");
    if (!oterm) oterm = inst->findITerm("Z");
    if (oterm && oterm->getNet()) {
      for (auto* iterm : oterm->getNet()->getITerms()) {
        if (iterm == oterm) continue;
        auto* driven_inst = iterm->getInst();
        std::string pin_name = iterm->getMTerm()->getName();
        if (pin_name == "CLK" || pin_name == "CK") {
          buf.driven_ffs.push_back(driven_inst->getName());
        }
      }
    }

    buf.current_delay_ps = getBufferDelay(buf.cell_name);

    int buf_idx = static_cast<int>(buffers_.size());
    for (const auto& ff : buf.driven_ffs) {
      ff_to_buf_[ff] = buf_idx;
    }
    buffers_.push_back(std::move(buf));
  }

  logger_->info(utl::CTS, CTS_BS_BUFS,
                "BufferSizing: {} clock buffers, {} driven FFs",
                buffers_.size(), ff_to_buf_.size());
}

void CtsBufferSizingLp::enumerateBufferSizes()
{
  // For each buffer, enumerate available sizing options from same tier
  // Build Pareto frontier: (delay, strength) non-dominated set

  // ASAP7 buffer library
  struct BufEntry { std::string suffix; int strength; };
  std::vector<BufEntry> asap7_bufs = {
    {"BUFx2",  2}, {"BUFx3",  3}, {"BUFx4",  4}, {"BUFx4f", 4},
    {"BUFx5",  5}, {"BUFx6f", 6}, {"BUFx8",  8}, {"BUFx10", 10},
    {"BUFx12", 12}, {"BUFx12f", 12}, {"BUFx16f", 16}, {"BUFx24", 24}
  };

  for (auto& buf : buffers_) {
    buf.size_options.clear();
    std::string tier_suffix = (buf.tier == "bottom") ? "_ASAP7_75t_R_bottom"
                                                      : "_ASAP7_75t_R_upper";
    // Deduplicate by strength (prefer non-f variant)
    std::map<int, BufSizeOption> by_strength;
    for (const auto& entry : asap7_bufs) {
      std::string master = entry.suffix + tier_suffix;
      float delay = getBufferDelay(master);
      float cap = entry.strength * 0.125f;  // fF estimate

      BufSizeOption opt{master, delay, entry.strength, cap};
      auto it = by_strength.find(entry.strength);
      if (it == by_strength.end()) {
        by_strength[entry.strength] = opt;
      } else {
        // Prefer non-f variant
        bool existing_f = it->second.master_name.find('f') != std::string::npos;
        bool new_f = master.find("BUFx") != std::string::npos
                  && master.substr(master.find("BUFx") + 4, 1) != "_"
                  && master.find('f') < master.find("_ASAP7");
        if (existing_f && !new_f) {
          by_strength[entry.strength] = opt;
        }
      }
    }

    for (const auto& [_, opt] : by_strength) {
      buf.size_options.push_back(opt);
    }
    // Sort by delay descending (slowest first)
    std::sort(buf.size_options.begin(), buf.size_options.end(),
              [](const BufSizeOption& a, const BufSizeOption& b) {
                return a.delay_ps > b.delay_ps;
              });
  }
}

std::vector<BufSizingChange> CtsBufferSizingLp::solve(float hold_weight,
                                                       float reg_weight)
{
  const int n_bufs = static_cast<int>(buffers_.size());
  if (n_bufs == 0) {
    return {};
  }

  // Create GLOP solver
  std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("GLOP"));
  if (!solver) {
    logger_->error(utl::CTS, CTS_BS_LP,
                   "BufferSizing: Failed to create GLOP solver");
    return {};
  }
  const double inf = solver->infinity();

  // Variables: d_k (delay per buffer)
  std::vector<MPVariable*> d_vars(n_bufs);
  std::vector<float> current_d(n_bufs);
  for (int k = 0; k < n_bufs; ++k) {
    current_d[k] = buffers_[k].current_delay_ps;
    float d_min = current_d[k], d_max = current_d[k];
    if (!buffers_[k].size_options.empty()) {
      d_min = buffers_[k].size_options.back().delay_ps;   // fastest (smallest)
      d_max = buffers_[k].size_options.front().delay_ps;  // slowest (largest)
    }
    d_vars[k] = solver->MakeNumVar(d_min, d_max, "d_" + std::to_string(k));
  }

  // Soft setup violation variables (one per edge)
  const int n_edges = static_cast<int>(edges_.size());
  std::vector<MPVariable*> v_vars(n_edges);
  for (int e = 0; e < n_edges; ++e) {
    v_vars[e] = solver->MakeNumVar(0.0, inf, "v_" + std::to_string(e));
  }

  // Objective: penalty for violations + regularization
  MPObjective* obj = solver->MutableObjective();
  obj->SetMinimization();

  // Penalty for setup violations
  for (int e = 0; e < n_edges; ++e) {
    obj->SetCoefficient(v_vars[e], 10.0);  // high penalty
  }

  // Setup/hold improvement weighted in objective via d_k coefficients
  for (int e = 0; e < n_edges; ++e) {
    const auto& edge = edges_[e];
    auto it_i = ff_to_buf_.find(edge.from_ff);
    auto it_j = ff_to_buf_.find(edge.to_ff);
    if (it_i == ff_to_buf_.end() || it_j == ff_to_buf_.end()) continue;
    int ki = it_i->second;
    int kj = it_j->second;
    if (ki == kj) continue;

    if (edge.slack_setup_ps < 0) {
      float w = std::abs(edge.slack_setup_ps);
      obj->SetCoefficient(d_vars[ki], obj->GetCoefficient(d_vars[ki]) + w);
      obj->SetCoefficient(d_vars[kj], obj->GetCoefficient(d_vars[kj]) - w);
    }
  }

  // Constraints
  for (int e = 0; e < n_edges; ++e) {
    const auto& edge = edges_[e];
    auto it_i = ff_to_buf_.find(edge.from_ff);
    auto it_j = ff_to_buf_.find(edge.to_ff);
    if (it_i == ff_to_buf_.end() || it_j == ff_to_buf_.end()) continue;
    int ki = it_i->second;
    int kj = it_j->second;
    if (ki == kj) continue;

    // Setup: d_i - d_j - v_e ≤ slack_setup + (current_d_i - current_d_j)
    {
      float rhs = edge.slack_setup_ps + current_d[ki] - current_d[kj];
      MPConstraint* ct = solver->MakeRowConstraint(-inf, rhs, "");
      ct->SetCoefficient(d_vars[ki], 1.0);
      ct->SetCoefficient(d_vars[kj], -1.0);
      ct->SetCoefficient(v_vars[e], -1.0);
    }

    // Hold: d_j - d_i ≤ slack_hold + (current_d_j - current_d_i)  (HARD)
    {
      float rhs = edge.slack_hold_ps + current_d[kj] - current_d[ki];
      MPConstraint* ct = solver->MakeRowConstraint(-inf, rhs, "");
      ct->SetCoefficient(d_vars[kj], 1.0);
      ct->SetCoefficient(d_vars[ki], -1.0);
    }
  }

  logger_->info(utl::CTS, CTS_BS_LP,
                "BufferSizing: {} buffers, {} edges, solving...",
                n_bufs, n_edges);

  // Solve
  auto t0 = std::chrono::high_resolution_clock::now();
  const MPSolver::ResultStatus status = solver->Solve();
  auto t1 = std::chrono::high_resolution_clock::now();
  double elapsed = std::chrono::duration<double>(t1 - t0).count();

  if (status != MPSolver::OPTIMAL && status != MPSolver::FEASIBLE) {
    logger_->warn(utl::CTS, CTS_BS_RESULT,
                  "BufferSizing: LP failed (status={})", static_cast<int>(status));
    return {};
  }

  // Map continuous delays to discrete buffer sizes
  std::vector<BufSizingChange> changes;
  for (int k = 0; k < n_bufs; ++k) {
    float opt_delay = static_cast<float>(d_vars[k]->solution_value());
    const auto& options = buffers_[k].size_options;
    if (options.empty()) continue;

    // Find closest discrete option
    const BufSizeOption* best = &options[0];
    float best_diff = std::abs(options[0].delay_ps - opt_delay);
    for (const auto& opt : options) {
      float diff = std::abs(opt.delay_ps - opt_delay);
      if (diff < best_diff) {
        best_diff = diff;
        best = &opt;
      }
    }

    if (best->master_name != buffers_[k].cell_name) {
      BufSizingChange change;
      change.inst_name = buffers_[k].inst_name;
      change.old_master = buffers_[k].cell_name;
      change.new_master = best->master_name;
      change.old_delay_ps = current_d[k];
      change.new_delay_ps = best->delay_ps;
      changes.push_back(change);
    }
  }

  logger_->info(utl::CTS, CTS_BS_RESULT,
                "BufferSizing: Solved in {:.2f}s, {} buffers resized",
                elapsed, changes.size());

  return changes;
}

void CtsBufferSizingLp::applyChanges(const std::vector<BufSizingChange>& changes)
{
  int applied = 0;
  int pin_remapped = 0;
  for (const auto& change : changes) {
    auto* inst = block_->findInst(change.inst_name.c_str());
    if (!inst) continue;
    auto* new_master = block_->getDb()->findMaster(change.new_master.c_str());
    if (!new_master) continue;

    // Pin-remap swap for heterogeneous 3D.
    // BUF_X4_upper (output pin Z) vs BUFx4_ASAP7_75t_R_upper (output pin Y).
    // swapMaster checks master mterm names — fails if Z != Y.
    // Fix: detect pin mismatch → destroy instance → recreate with new master
    // → reconnect nets by IO type (input A→A, output Z→Y).
    auto* old_master = inst->getMaster();

    // Try normal swap first (fast path, same-pin-name cells)
    if (inst->swapMaster(new_master)) {
      ++applied;
      continue;
    }

    // swapMaster failed (likely pin name mismatch). Do manual remap.
    // 1. Save net connections by IO type
    odb::dbNet* in_net = nullptr;
    odb::dbNet* out_net = nullptr;
    for (auto* iterm : inst->getITerms()) {
      auto* mt = iterm->getMTerm();
      if (mt->getSigType() != odb::dbSigType::SIGNAL) continue;
      if (mt->getIoType() == odb::dbIoType::INPUT && !in_net) {
        in_net = iterm->getNet();
      } else if (mt->getIoType() == odb::dbIoType::OUTPUT && !out_net) {
        out_net = iterm->getNet();
      }
    }

    if (!in_net && !out_net) continue;  // no signal connections, skip

    // 2. Save placement
    int x, y;
    inst->getLocation(x, y);
    auto orient = inst->getOrient();
    auto status = inst->getPlacementStatus();
    std::string inst_name = inst->getName();

    // 3. Destroy old instance
    odb::dbInst::destroy(inst);

    // 4. Create new instance with new master
    auto* new_inst = odb::dbInst::create(block_, new_master, inst_name.c_str());
    if (!new_inst) {
      logger_->warn(utl::CTS, CTS_BS_RESULT + 2,
          "Pin-remap: failed to recreate {}", inst_name);
      continue;
    }

    // 5. Restore placement
    new_inst->setLocation(x, y);
    new_inst->setOrient(orient);
    new_inst->setPlacementStatus(status);

    // 6. Connect nets by IO type (input→input, output→output)
    for (auto* iterm : new_inst->getITerms()) {
      auto* mt = iterm->getMTerm();
      if (mt->getSigType() != odb::dbSigType::SIGNAL) continue;
      if (mt->getIoType() == odb::dbIoType::INPUT && in_net) {
        iterm->connect(in_net);
      } else if (mt->getIoType() == odb::dbIoType::OUTPUT && out_net) {
        iterm->connect(out_net);
      }
    }

    ++applied;
    ++pin_remapped;
  }
  logger_->info(utl::CTS, CTS_BS_RESULT + 1,
                "BufferSizing: Applied {} master swaps ({} pin-remapped)",
                applied, pin_remapped);
}

}  // namespace cts
