// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "cts/TritonCTS.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <ranges>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Clock.h"
#include "Cts3DDatabase.h"        // Added for 3D tier database
#include "CtsSkewLpSolver.h"     // C++ LP-TNS solver
#include "BufSizingLpSolver.h"   // C++ buffer sizing LP
// CtsBufferSizingLp.h removed — replaced by BufSizingLpSolver (V53_FM_BUF)
#include "ClockLatencyEstimator.h"  // per-FF physical achievability
#include "CtsOptions.h"
#include "HTreeBuilder.h"
#include "LatencyBalancer.h"
#include "TechChar.h"
#include "TreeBuilder.h"
#include "VerilogFFExtractor.h"
#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "est/EstimateParasitics.h"
#include "odb/db.h"
#include "odb/dbSet.h"
#include "odb/dbShape.h"
#include "odb/dbTypes.h"
#include "odb/geom.h"
#include "rsz/Resizer.hh"
#include "sta/Clock.hh"
#include "sta/Fuzzy.hh"
#include "sta/Graph.hh"
#include "sta/GraphDelayCalc.hh"
#include "sta/Liberty.hh"
#include "sta/Network.hh"
#include "sta/NetworkClass.hh"
#include "sta/PathAnalysisPt.hh"
#include "sta/PathEnd.hh"
#include "sta/PatternMatch.hh"
#include "sta/Sdc.hh"
#include "sta/Vector.hh"
#include "stt/SteinerTreeBuilder.h"
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

TritonCTS::TritonCTS(utl::Logger* logger,
                     odb::dbDatabase* db,
                     sta::dbNetwork* network,
                     sta::dbSta* sta,
                     stt::SteinerTreeBuilder* st_builder,
                     rsz::Resizer* resizer,
                     est::EstimateParasitics* estimate_parasitics)
{
  logger_ = logger;
  db_ = db;
  network_ = network;
  openSta_ = sta;
  resizer_ = resizer;
  estimate_parasitics_ = estimate_parasitics;

  options_ = new CtsOptions(logger_, st_builder);
}

TritonCTS::~TritonCTS()
{
  delete options_;
}

void TritonCTS::runTritonCts()
{
  odb::dbChip* chip = db_->getChip();
  odb::dbBlock* block = chip->getBlock();
  options_->addOwner(block);

  setupCharacterization();

  // Initialize 3D tier database from ODB
  // Preserve pre-loaded cts3dDb_ (e.g., with skew targets)
  const bool hadPreloadedDb = (cts3dDb_ != nullptr);
  if (!hadPreloadedDb) {
    cts3dDb_ = std::make_unique<Cts3DDatabase>(db_, logger_);
  }
  cts3dDb_->populate();

  // Copy wire RC from TechChar to 3D database
  // Initially same for both tiers, can be split later if needed
  const double wire_res = techChar_->getResPerDBU();
  const double wire_cap = techChar_->getCapPerDBU();
  cts3dDb_->setWireRC(0, wire_res, wire_cap);  // bottom tier
  cts3dDb_->setWireRC(1, wire_res, wire_cap);  // upper tier

  // Auto-detect bottom/upper buffer pair from buf_list.
  // When buf_list contains cells from different technology nodes (e.g., NG45
  // "BUF_X4_bottom" and ASAP7 "BUF_X4_upper"), suffix substitution alone in
  // getBufferForTier() fails.  Register an explicit pair so that the correct
  // cell is selected for each tier during tree building.
  {
    // Lambda to check string suffix (hasSuffix is private in Cts3DDatabase)
    auto hasSfx = [](const std::string& s, const std::string& sfx) {
      return s.size() >= sfx.size()
             && s.compare(s.size() - sfx.size(), sfx.size(), sfx) == 0;
    };
    std::string bottomBuf, upperBuf;
    for (const auto& buf : options_->getBufferList()) {
      if (hasSfx(buf, "_bottom")) {
        bottomBuf = buf;
      } else if (hasSfx(buf, "_upper")) {
        upperBuf = buf;
      }
    }
    if (!bottomBuf.empty() && !upperBuf.empty()) {
      cts3dDb_->setTierBufferPair(bottomBuf, upperBuf);
    }
  }

  findClockRoots();
  populateTritonCTS();
  if (builders_.empty()) {
    logger_->warn(CTS, 82, "No valid clock nets in the design.");
  } else {
    checkCharacterization();
    buildClockTrees();
    writeDataToDb();
    setAllClocksPropagated();
    if (options_->getRepairClockNets()) {
      repairClockNets();
    }
    balanceMacroRegisterLatencies();
  }

  // reset
  techChar_.reset();
  // Only reset 3D DB if no pre-loaded skew targets
  if (!cts3dDb_ || !cts3dDb_->hasSkewTargets()) {
    cts3dDb_.reset();  // Cleanup 3D tier database
  }
  builders_.clear();
  staClockNets_.clear();
  visitedClockNets_.clear();
  inst2clkbuf_.clear();
  driver2subnet_.clear();
  numberOfClocks_ = 0;
  numClkNets_ = 0;
  numFixedNets_ = 0;
  dummyLoadIndex_ = 0;
  rootBuffers_.clear();
  sinkBuffers_.clear();
  regTreeRootBufIndex_ = 0;
  delayBufIndex_ = 0;
  options_->removeOwner();
}

TreeBuilder* TritonCTS::addBuilder(CtsOptions* options,
                                   Clock& net,
                                   odb::dbNet* topInputNet,
                                   TreeBuilder* parent,
                                   utl::Logger* logger,
                                   odb::dbDatabase* db)
{
  auto builder
      = std::make_unique<HTreeBuilder>(options, net, parent, logger, db);

  builder->setTopInputNet(topInputNet);
  builders_.emplace_back(std::move(builder));
  return builders_.back().get();
}

int TritonCTS::getBufferFanoutLimit(const std::string& bufferName)
{
  int fanout = std::numeric_limits<int>::max();
  float tempFanout;
  bool existMaxFanout;

  // Check if top instance has fanout limit
  sta::Cell* top_cell = network_->cell(network_->topInstance());
  openSta_->sdc()->fanoutLimit(
      top_cell, sta::MinMax::max(), tempFanout, existMaxFanout);
  if (existMaxFanout) {
    fanout = std::min(fanout, (int) tempFanout);
  }

  odb::dbMaster* bufferMaster = db_->findMaster(bufferName.c_str());
  sta::Cell* bufferCell = network_->dbToSta(bufferMaster);
  sta::Port* buffer_port = nullptr;
  for (odb::dbMTerm* mterm : bufferMaster->getMTerms()) {
    odb::dbSigType sig_type = mterm->getSigType();
    if (sig_type == odb::dbSigType::GROUND
        || sig_type == odb::dbSigType::POWER) {
      continue;
    }
    odb::dbIoType io_type = mterm->getIoType();
    if (io_type == odb::dbIoType::OUTPUT) {
      buffer_port = network_->dbToSta(mterm);
      break;
    }
  }
  if (buffer_port == nullptr) {
    return (existMaxFanout) ? fanout : 0;
  }

  openSta_->sdc()->fanoutLimit(
      buffer_port, sta::MinMax::max(), tempFanout, existMaxFanout);
  if (existMaxFanout) {
    fanout = std::min(fanout, (int) tempFanout);
  }

  openSta_->sdc()->fanoutLimit(
      bufferCell, sta::MinMax::max(), tempFanout, existMaxFanout);
  if (existMaxFanout) {
    fanout = std::min(fanout, (int) tempFanout);
  }

  sta::LibertyPort* port = network_->libertyPort(buffer_port);
  port->fanoutLimit(sta::MinMax::max(), tempFanout, existMaxFanout);
  if (existMaxFanout) {
    fanout = std::min(fanout, (int) tempFanout);
  } else {
    port->libertyLibrary()->defaultMaxFanout(tempFanout, existMaxFanout);
    if ((existMaxFanout)) {
      fanout = std::min(fanout, (int) tempFanout);
    }
  }
  return fanout == std::numeric_limits<int>::max() ? 0 : fanout;
}

void TritonCTS::setupCharacterization()
{
  // Check if CTS library is valid
  if (options_->isCtsLibrarySet()) {
    sta::Library* lib = network_->findLibrary(options_->getCtsLibrary());
    if (lib == nullptr) {
      logger_->error(CTS,
                     209,
                     "Library {} cannot be found because it is not "
                     "loaded or name is incorrect",
                     options_->getCtsLibrary());
    } else {
      logger_->info(CTS,
                    210,
                    "Clock buffers will be chosen from library {}",
                    options_->getCtsLibrary());
    }
  }

  block_ = db_->getChip()->getBlock();
  options_->setDbUnits(block_->getDbUnitsPerMicron());

  openSta_->checkFanoutLimitPreamble();
  // Finalize root/sink buffers
  std::string rootBuffer = selectRootBuffer(rootBuffers_);
  options_->setRootBuffer(rootBuffer);
  std::string sinkBuffer = selectSinkBuffer(sinkBuffers_);
  options_->setSinkBuffer(sinkBuffer);

  int sinkMaxFanout = getBufferFanoutLimit(sinkBuffer);
  int rootMaxFanout = getBufferFanoutLimit(rootBuffer);

  if (rootMaxFanout && (options_->getNumMaxLeafSinks() > rootMaxFanout)) {
    options_->setMaxFanout(rootMaxFanout);
  }

  if (sinkMaxFanout) {
    options_->limitSinkClusteringSizes(sinkMaxFanout);
    if (sinkMaxFanout < options_->getMaxFanout()) {
      options_->setMaxFanout(sinkMaxFanout);
    }
  }

  double maxWlMicrons
      = resizer_->findMaxWireLength(/* don't issue error */ false) * 1e+6;
  options_->setMaxWl(block_->micronsToDbu(maxWlMicrons));

  // A new characteriztion is always created.
  techChar_ = std::make_unique<TechChar>(options_,
                                         db_,
                                         openSta_,
                                         resizer_,
                                         estimate_parasitics_,
                                         network_,
                                         logger_);
  techChar_->create();

  // Also resets metrics everytime the setup is done
  options_->setNumSinks(0);
  options_->setNumBuffersInserted(0);
  options_->setNumClockRoots(0);
  options_->setNumClockSubnets(0);
}

void TritonCTS::checkCharacterization()
{
  std::unordered_set<std::string> visitedMasters;
  techChar_->forEachWireSegment([&](unsigned idx, const WireSegment& wireSeg) {
    for (int buf = 0; buf < wireSeg.getNumBuffers(); ++buf) {
      const std::string& master = wireSeg.getBufferMaster(buf);
      if (!visitedMasters.contains(master)) {
        if (masterExists(master)) {
          visitedMasters.insert(master);
        } else {
          logger_->error(CTS, 81, "Buffer {} is not in the loaded DB.", master);
        }
      }
    }
  });

  logger_->info(CTS,
                97,
                "Characterization used {} buffer(s) types.",
                visitedMasters.size());
}

void TritonCTS::findClockRoots()
{
  if (!options_->getClockNets().empty()) {
    logger_->info(CTS,
                  1,
                  "Running TritonCTS with user-specified clock roots: {}.",
                  options_->getClockNets());
  }
}

void TritonCTS::buildClockTrees()
{
  for (auto& builder : builders_) {
    builder->setTechChar(*techChar_);
    // Inject 3D tier database into each builder
    if (cts3dDb_) {
      builder->setCts3DDatabase(*cts3dDb_);
    }
    builder->setDb(db_);
    builder->setLogger(logger_);
    builder->initBlockages();
    builder->run();
  }
}

void TritonCTS::initOneClockTree(odb::dbNet* driverNet,
                                 odb::dbNet* clkInputNet,
                                 const std::string& sdcClockName,
                                 TreeBuilder* parent)
{
  TreeBuilder* clockBuilder = nullptr;
  std::vector<odb::dbNet*> skipNets = options_->getSkipNets();
  if (driverNet->isSpecial()) {
    logger_->info(
        CTS, 116, "Special net \"{}\" skipped.", driverNet->getName());
  } else if (std::ranges::find(skipNets, driverNet) != skipNets.end()) {
    logger_->warn(CTS,
                  44,
                  "Skipping net {}, specified by the user...",
                  driverNet->getName());
  } else {
    clockBuilder = initClock(driverNet, clkInputNet, sdcClockName, parent);
  }
  if (clockBuilder != nullptr && net2builder_[clkInputNet] == nullptr) {
    net2builder_[clkInputNet] = clockBuilder;
  }
  // Treat gated clocks as separate clock trees
  // TODO: include sinks from gated clocks together with other sinks and build
  // one clock tree
  visitedClockNets_.insert(driverNet);
  odb::dbITerm* driver = driverNet->getFirstOutput();
  odb::dbSet<odb::dbITerm> iterms = driverNet->getITerms();
  for (odb::dbITerm* iterm : iterms) {
    if (iterm != driver && iterm->isInputSignal()) {
      if (!isSink(iterm)) {
        odb::dbITerm* outputPin = getSingleOutput(iterm->getInst(), iterm);
        if (outputPin && outputPin->getNet()) {
          odb::dbNet* outputNet = outputPin->getNet();
          if (visitedClockNets_.find(outputNet) == visitedClockNets_.end()
              && !openSta_->sdc()->isLeafPinClock(
                  network_->dbToSta(outputPin))) {
            if (clockBuilder == nullptr
                && net2builder_[clkInputNet] != nullptr) {
              initOneClockTree(outputNet,
                               clkInputNet,
                               sdcClockName,
                               net2builder_[clkInputNet]);
            } else {
              initOneClockTree(
                  outputNet, clkInputNet, sdcClockName, clockBuilder);
            }
          }
        }
      }
    }
  }
}

// Iterative version of countSinksPostDbWrite.
// Original recursive version caused stack overflow on large designs
// (ariane133: 5,800 nets level 9, bp_quad: 50,000+ nets level 12+).
// Converted to explicit std::stack on heap — identical results, no depth limit.
void TritonCTS::countSinksPostDbWrite(
    TreeBuilder* builder,
    odb::dbNet* net,
    unsigned& sinks_cnt,
    unsigned& leafSinks,
    unsigned currWireLength,
    double& sinkWireLength,
    int& minDepth,
    int& maxDepth,
    int depth,
    bool fullTree,
    const std::unordered_set<odb::dbITerm*>& sinks,
    const std::unordered_set<odb::dbInst*>& dummies)
{
  struct StackFrame {
    odb::dbNet* net;
    unsigned wireLength;
    int depth;
  };

  std::stack<StackFrame> work;
  work.push({net, currWireLength, depth});

  while (!work.empty()) {
    auto [curNet, curWireLen, curDepth] = work.top();
    work.pop();

    // Find driver position
    int driverX = 0, driverY = 0;
    for (odb::dbITerm* iterm : curNet->getITerms()) {
      if (iterm->getIoType() != odb::dbIoType::INPUT) {
        iterm->getAvgXY(&driverX, &driverY);
        break;
      }
    }
    for (odb::dbBTerm* bterm : curNet->getBTerms()) {
      if (bterm->getIoType() == odb::dbIoType::INPUT) {
        for (odb::dbBPin* pin : bterm->getBPins()) {
          odb::dbPlacementStatus status = pin->getPlacementStatus();
          if (status == odb::dbPlacementStatus::NONE
              || status == odb::dbPlacementStatus::UNPLACED) {
            continue;
          }
          for (odb::dbBox* box : pin->getBoxes()) {
            if (box) {
              driverX = box->xMin();
              driverY = box->yMin();
              break;
            }
          }
          break;
        }
      }
    }

    // Process input iterms
    for (odb::dbITerm* iterm : curNet->getITerms()) {
      if (iterm->getIoType() != odb::dbIoType::INPUT) continue;

      int receiverX, receiverY;
      iterm->getAvgXY(&receiverX, &receiverY);
      unsigned dist = abs(driverX - receiverX) + abs(driverY - receiverY);
      odb::dbInst* inst = iterm->getInst();

      bool terminate = fullTree
                           ? (sinks.find(iterm) != sinks.end())
                           : !builder->isAnyTreeBuffer(getClockFromInst(inst));
      odb::dbITerm* outputPin = inst->getFirstOutput();
      bool trueSink = true;

      if (outputPin && outputPin->getNet() == curNet) {
        terminate = true;
        trueSink = false;
      }

      if (!terminate && inst) {
        if (isMacroBlockInst(inst)) {
          terminate = true;
          trueSink = false;
        } else {
          sta::Cell* masterCell = network_->dbToSta(inst->getMaster());
          if (masterCell) {
            sta::LibertyCell* libCell = network_->libertyCell(masterCell);
            if (libCell && libCell->hasSequentials()) {
              terminate = true;
              trueSink = false;
            }
          }
        }
      }

      if (!terminate) {
        if (outputPin && outputPin->getNet() != nullptr) {
          // Push downstream net onto stack (was: recursive call)
          work.push({outputPin->getNet(), curWireLen + dist, curDepth + 1});
        } else {
          std::string cellType = "Complex cell";
          sta::Cell* masterCell = network_->dbToSta(inst->getMaster());
          if (masterCell) {
            sta::LibertyCell* libCell = network_->libertyCell(masterCell);
            if (libCell) {
              if (libCell->isInverter()) cellType = "Inverter";
              else if (libCell->isBuffer()) cellType = "Buffer";
            }
          }
          if (dummies.find(inst) == dummies.end()) {
            logger_->info(CTS, 121, "{} '{}' has unconnected output pin.",
                          cellType, inst->getName());
          }
        }
        if (builder->isLeafBuffer(getClockFromInst(inst))) {
          leafSinks++;
        }
      } else if (trueSink) {
        sinks_cnt++;
        double currSinkWl
            = (dist + curWireLen) / double(options_->getDbUnits());
        sinkWireLength += currSinkWl;
        maxDepth = std::max(curDepth, maxDepth);
        if ((minDepth > 0 && curDepth < minDepth) || (minDepth == 0)) {
          minDepth = curDepth;
        }
      }
    }
  }
}

ClockInst* TritonCTS::getClockFromInst(odb::dbInst* inst)
{
  auto it = inst2clkbuf_.find(inst);
  return it != inst2clkbuf_.end() ? it->second : nullptr;
}

void TritonCTS::writeDataToDb()
{
  std::set<odb::dbNet*> clkLeafNets;
  std::unordered_set<odb::dbInst*> clkDummies;

  for (auto& builder : builders_) {
    writeClockNetsToDb(builder.get(), clkLeafNets);
    if (options_->getApplyNdr() != CtsOptions::NdrStrategy::NONE) {
      writeClockNDRsToDb(builder.get());
    }
    if (options_->dummyLoadEnabled()) {
      int nDummies = writeDummyLoadsToDb(builder->getClock(), clkDummies);
      builder->setNDummies(nDummies);
    }
  }

  for (auto& builder : builders_) {
    odb::dbNet* topClockNet = builder->getClock().getNetObj();
    unsigned sinkCount = 0;
    unsigned leafSinks = 0;
    double allSinkDistance = 0.0;
    int minDepth = 0;
    int maxDepth = 0;
    bool reportFullTree
        = !builder->getParent() && !builder->getChildren().empty();

    std::unordered_set<odb::dbITerm*> sinks;
    builder->getClock().forEachSink([&sinks](const ClockInst& inst) {
      sinks.insert(inst.getDbInputPin());
    });
    if (sinks.size() < 2) {
      logger_->info(
          CTS, 124, "Clock net \"{}\"", builder->getClock().getName());
      logger_->info(CTS, 125, " Sinks {}", sinks.size());
    } else {
      countSinksPostDbWrite(builder.get(),
                            topClockNet,
                            sinkCount,
                            leafSinks,
                            0,
                            allSinkDistance,
                            minDepth,
                            maxDepth,
                            0,
                            reportFullTree,
                            sinks,
                            clkDummies);
      logger_->info(CTS, 98, "Clock net \"{}\"", builder->getClock().getName());
      logger_->info(CTS, 99, " Sinks {}", sinkCount);
      logger_->info(CTS, 100, " Leaf buffers {}", leafSinks);
      if (sinkCount > 0) {
        double avgWL = allSinkDistance / sinkCount;
        logger_->info(CTS, 101, " Average sink wire length {:.2f} um", avgWL);
      }
      logger_->info(CTS, 102, " Path depth {} - {}", minDepth, maxDepth);
      if (options_->dummyLoadEnabled()) {
        logger_->info(
            CTS, 207, " Dummy loads inserted {}", builder->getNDummies());
      }
    }
  }
}

void TritonCTS::forEachBuilder(
    const std::function<void(const TreeBuilder*)>& func) const
{
  for (const auto& builder : builders_) {
    func(builder.get());
  }
}

void TritonCTS::reportCtsMetrics()
{
  std::string filename = options_->getMetricsFile();

  if (!filename.empty()) {
    std::ofstream file(filename.c_str());

    if (!file.is_open()) {
      logger_->error(
          CTS, 87, "Could not open output metric file {}.", filename.c_str());
    }

    file << "Total number of Clock Roots: " << options_->getNumClockRoots()
         << ".\n";
    file << "Total number of Buffers Inserted: "
         << options_->getNumBuffersInserted() << ".\n";
    file << "Total number of Clock Subnets: " << options_->getNumClockSubnets()
         << ".\n";
    file << "Total number of Sinks: " << options_->getNumSinks() << ".\n";

    file << "Buffers used:\n";
    for (const auto& [master, count] : options_->getBufferCount()) {
      file << "  " << master->getName() << ": " << count << "\n";
    }
    if (!options_->getDummyCount().empty()) {
      file << "Dummys used:\n";
      for (const auto& [master, count] : options_->getDummyCount()) {
        file << "  " << master->getName() << ": " << count << "\n";
      }
    }
    file.close();

  } else {
    logger_->report("Total number of Clock Roots: {}.",
                    options_->getNumClockRoots());
    logger_->report("Total number of Buffers Inserted: {}.",
                    options_->getNumBuffersInserted());
    logger_->report("Total number of Clock Subnets: {}.",
                    options_->getNumClockSubnets());
    logger_->report("Total number of Sinks: {}.", options_->getNumSinks());

    logger_->report("Cells used:");
    for (const auto& [master, count] : options_->getBufferCount()) {
      logger_->report("  {}: {}", master->getName(), count);
    }
    if (!options_->getDummyCount().empty()) {
      logger_->report("Dummys used:");
      for (const auto& [master, count] : options_->getDummyCount()) {
        logger_->report("  {}: {}", master->getName(), count);
      }
    }
  }
}

int TritonCTS::setClockNets(const char* names)
{
  odb::dbChip* chip = db_->getChip();
  odb::dbBlock* block = chip->getBlock();

  options_->setClockNets(names);
  std::stringstream ss(names);
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> nets(begin, end);

  std::vector<odb::dbNet*> netObjects;

  for (const std::string& name : nets) {
    odb::dbNet* net = block->findNet(name.c_str());
    bool netFound = false;
    if (net != nullptr) {
      // Since a set is unique, only the nets not found by dbSta are added.
      netObjects.push_back(net);
      netFound = true;
    } else {
      // User input was a pin, transform it into an iterm if possible
      odb::dbITerm* iterm = block->findITerm(name.c_str());
      if (iterm != nullptr) {
        net = iterm->getNet();
        if (net != nullptr) {
          // Since a set is unique, only the nets not found by dbSta are added.
          netObjects.push_back(net);
          netFound = true;
        }
      }
    }
    if (!netFound) {
      return 1;
    }
  }
  options_->setClockNetsObjs(netObjects);
  return 0;
}

void TritonCTS::setBufferList(const char* buffers)
{
  // Put the buffer list into a string vector
  std::stringstream ss(buffers);
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> bufferList(begin, end);
  // If the vector is empty, then the buffers are inferred
  if (bufferList.empty()) {
    const char* lib_name
        = options_->isCtsLibrarySet() ? options_->getCtsLibrary() : nullptr;
    resizer_->inferClockBufferList(lib_name, bufferList);
    options_->setBufferListInferred(true);
  } else {
    // Iterate the user-defined buffer list
    sta::Vector<sta::LibertyCell*> selected_buffers;
    for (const std::string& buffer : bufferList) {
      odb::dbMaster* buffer_master = db_->findMaster(buffer.c_str());
      if (buffer_master == nullptr) {
        logger_->error(
            CTS, 126, "No physical master cell found for buffer {}.", buffer);
      } else {
        // Get the buffer and add to the vector
        sta::Cell* master_cell = network_->dbToSta(buffer_master);
        if (master_cell) {
          sta::LibertyCell* lib_cell = network_->libertyCell(master_cell);
          selected_buffers.push_back(lib_cell);
        }
      }
    }
    // Add found buffer to RSZ
    resizer_->setClockBuffersList(selected_buffers);
  }
  options_->setBufferList(bufferList);
}

std::string TritonCTS::getRootBufferToString()
{
  std::ostringstream buffer_names;
  for (const auto& buf : rootBuffers_) {
    buffer_names << buf << " ";
  }
  return buffer_names.str();
}

void TritonCTS::setRootBuffer(const char* buffers)
{
  std::stringstream ss(buffers);
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> bufferList(begin, end);
  for (const std::string& buffer : bufferList) {
    if (db_->findMaster(buffer.c_str()) == nullptr) {
      logger_->error(
          CTS, 127, "No physical master cell found for buffer {}.", buffer);
    }
  }
  rootBuffers_ = std::move(bufferList);
}

std::string TritonCTS::selectRootBuffer(std::vector<std::string>& buffers)
{
  // if -root_buf is not specified, choose from the buffer list
  if (buffers.empty()) {
    buffers = options_->getBufferList();
  }

  if (buffers.size() == 1) {
    return buffers.front();
  }

  options_->setRootBufferInferred(true);
  // estimate wire cap for root buffer
  // assume sink buffer needs to drive clk buffers at two far ends of chip
  // at midpoint
  //
  //  --------------
  //  |      .     |
  //  |   ===x===  |
  //  |      .     |
  //  --------------
  odb::dbBlock* block = db_->getChip()->getBlock();
  odb::Rect coreArea = block->getCoreArea();
  float sinkWireLength
      = static_cast<float>(std::max(coreArea.dx(), coreArea.dy()))
        / block->getDbUnitsPerMicron();
  sta::Corner* corner = openSta_->cmdCorner();
  float rootWireCap = estimate_parasitics_->wireSignalCapacitance(corner) * 1e-6
                      * sinkWireLength / 2.0;
  std::string rootBuf = selectBestMaxCapBuffer(buffers, rootWireCap);
  return rootBuf;
}

void TritonCTS::setSinkBuffer(const char* buffers)
{
  std::stringstream ss(buffers);
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> bufferList(begin, end);
  sinkBuffers_ = std::move(bufferList);
}

std::string TritonCTS::selectSinkBuffer(std::vector<std::string>& buffers)
{
  // if -sink_clustering_buf is not specified, choose from the buffer list
  if (buffers.empty()) {
    buffers = options_->getBufferList();
  }

  if (buffers.size() == 1) {
    return buffers.front();
  }

  options_->setSinkBufferInferred(true);
  // estimate wire cap for sink buffer
  // assume sink buffer needs to drive clk buffers at two far ends of chip
  // to account for unknown pin caps
  //
  //  --------------
  //  |======x=====|
  //  |      .     |
  //  |----- .-----|
  //  |      .     |
  //  --------------
  odb::dbBlock* block = db_->getChip()->getBlock();
  odb::Rect coreArea = block->getCoreArea();
  float sinkWireLength
      = static_cast<float>(std::max(coreArea.dx(), coreArea.dy()))
        / block->getDbUnitsPerMicron();
  sta::Corner* corner = openSta_->cmdCorner();
  float sinkWireCap = estimate_parasitics_->wireSignalCapacitance(corner) * 1e-6
                      * sinkWireLength;

  std::string sinkBuf = selectBestMaxCapBuffer(buffers, sinkWireCap);
  // clang-format off
  debugPrint(logger_, CTS, "buffering", 1, "{} has been selected as sink "
             "buffer to drive sink wire cap of {:0.2e}", sinkBuf, sinkWireCap);
  // clang-format on
  return sinkBuf;
}

// pick the smallest buffer that can drive total cap
// if no such buffer exists, pick one that has the largest max cap
std::string TritonCTS::selectBestMaxCapBuffer(
    const std::vector<std::string>& buffers,
    float totalCap)
{
  std::string bestBuf, nextBestBuf;
  float bestArea = std::numeric_limits<float>::max();
  float bestCap = 0.0;

  for (const std::string& name : buffers) {
    odb::dbMaster* master = db_->findMaster(name.c_str());
    if (master == nullptr) {
      logger_->error(
          CTS, 117, "Physical master could not be found for cell '{}'", name);
    }
    sta::Cell* masterCell = network_->dbToSta(master);
    sta::LibertyCell* libCell = network_->libertyCell(masterCell);
    if (libCell == nullptr) {
      logger_->error(
          CTS, 112, "Liberty cell could not be found for cell '{}'", name);
    }
    sta::LibertyPort *in, *out;
    libCell->bufferPorts(in, out);
    float area = libCell->area();
    float maxCap = 0.0;
    bool maxCapExists = false;
    out->capacitanceLimit(sta::MinMax::max(), maxCap, maxCapExists);
    // clang-format off
    debugPrint(logger_, CTS, "buffering", 1, "{} has cap limit:{}"
               " vs. total cap:{}, derate:{}", name,
               maxCap * (float) options_->getSinkBufferMaxCapDerate(), totalCap,
               options_->getSinkBufferMaxCapDerate());
    // clang-format on
    if (maxCapExists
        && ((maxCap * options_->getSinkBufferMaxCapDerate()) > totalCap)
        && area < bestArea) {
      bestBuf = name;
      bestArea = area;
    }
    if (maxCap > bestCap) {
      nextBestBuf = name;
      bestCap = maxCap;
    }
  }

  if (bestBuf.empty()) {
    bestBuf = std::move(nextBestBuf);
  }

  return bestBuf;
}

// db functions

void TritonCTS::cloneClockGaters(odb::dbNet* clkNet)
{
  odb::dbITerm* driver = clkNet->getFirstOutput();
  std::vector<int> xs;
  std::vector<int> ys;
  std::map<odb::Point, std::vector<odb::dbITerm*>> point2pin;
  odb::dbSet<odb::dbITerm> iterms = clkNet->getITerms();
  for (odb::dbITerm* iterm : iterms) {
    if (iterm != driver && iterm->isInputSignal()) {
      int TestX, TestY;
      iterm->getAvgXY(&TestX, &TestY);
      xs.push_back(TestX);
      ys.push_back(TestY);
      point2pin[{TestX, TestY}].push_back(iterm);
      if (isSink(iterm)) {
        continue;
      }
      odb::dbITerm* outputPin = getSingleOutput(iterm->getInst(), iterm);
      if (!outputPin || !outputPin->getNet()) {
        continue;
      }
      odb::dbInst* icg = iterm->getInst();
      odb::dbNet* outputNet = outputPin->getNet();
      sta::Cell* masterCell = network_->dbToSta(icg->getMaster());
      sta::LibertyCell* libertyCell = network_->libertyCell(masterCell);

      if (!libertyCell) {
        continue;
      }
      // Clock tree buffers or inverters
      if (libertyCell->isInverter() || libertyCell->isBuffer()) {
        continue;
      }
      cloneClockGaters(outputNet);
    }
  }
  if (!driver) {
    return;
  }

  // xs is empty if the fanout is a bterm
  if (isSink(driver) || driver->getInst()->isFixed()
      || driver->getInst()->isPad() || xs.empty()) {
    return;
  }

  int drvrX, drvrY;
  driver->getAvgXY(&drvrX, &drvrY);
  point2pin[{drvrX, drvrY}].push_back(driver);
  stt::Tree ftree
      = options_->getSttBuilder()->makeSteinerTree(clkNet, xs, ys, 0);
  findLongEdges(ftree, {drvrX, drvrY}, point2pin);
}

void TritonCTS::findLongEdges(
    stt::Tree& clkSteiner,
    odb::Point driverPt,
    std::map<odb::Point, std::vector<odb::dbITerm*>>& point2pin)
{
  const int threshold = options_->getMaxWl();
  debugPrint(
      logger_, CTS, "clock gate cloning", 1, "Threshold = {}", threshold);

  std::map<int, int> iterm2cluster;
  std::vector<std::vector<int>> clusters;
  odb::dbNet* icgNet = point2pin[driverPt][0]->getNet();
  odb::dbITerm* icgTerm = icgNet->getFirstOutput();
  std::string icgName = icgTerm->getInst()->getName();

  for (int b = 0; b < clkSteiner.branchCount(); b++) {
    const stt::Branch branch = clkSteiner.branch[b];
    const stt::Branch* neighbor = &clkSteiner.branch[branch.n];
    odb::Point branchPt = {branch.x, branch.y};

    odb::Point neighborPt = {neighbor->x, neighbor->y};
    int64_t dist = odb::Point::manhattanDistance(branchPt, neighborPt);
    const int clusterFrom
        = iterm2cluster.find(b) == iterm2cluster.end() ? -1 : iterm2cluster[b];
    const int clusterTo = iterm2cluster.find(branch.n) == iterm2cluster.end()
                              ? -1
                              : iterm2cluster[branch.n];

    if (b == branch.n) {
      continue;
    }

    if (dist >= threshold) {
      if (clusterFrom == -1) {
        int newClusterID = clusters.size();
        iterm2cluster[b] = newClusterID;
        clusters.push_back({b});
      }
      if (clusterTo == -1) {
        int newClusterID = clusters.size();
        iterm2cluster[branch.n] = newClusterID;
        clusters.push_back({branch.n});
      }
      continue;
    }

    if (clusterFrom != -1 && clusterTo != -1) {
      int mantainedCLuster
          = (clusters[clusterFrom].size() >= clusters[clusterTo].size())
                ? clusterFrom
                : clusterTo;
      int removedCLuster
          = (clusters[clusterFrom].size() < clusters[clusterTo].size())
                ? clusterFrom
                : clusterTo;

      clusters[mantainedCLuster].insert(clusters[mantainedCLuster].end(),
                                        clusters[removedCLuster].begin(),
                                        clusters[removedCLuster].end());
      for (int point : clusters[removedCLuster]) {
        iterm2cluster[point] = mantainedCLuster;
      }
      clusters[removedCLuster].clear();

    } else if (clusterFrom != -1) {
      iterm2cluster[branch.n] = clusterFrom;
      clusters[clusterFrom].push_back(branch.n);
    } else if (clusterTo != -1) {
      iterm2cluster[b] = clusterTo;
      clusters[clusterTo].push_back(b);
    } else {
      int newClusterID = clusters.size();
      iterm2cluster[b] = newClusterID;
      iterm2cluster[branch.n] = newClusterID;
      clusters.push_back({b, branch.n});
    }
  }

  // Find closest cluster to original ICG
  int driverClusterID = -1;
  int64_t minDist2Driver = std::numeric_limits<int64_t>::max();
  int validClusters = 0;
  for (int n = 0; n < clusters.size(); n++) {
    const std::vector<int>& cluster = clusters[n];
    if (cluster.empty()) {
      continue;
    }
    bool validCluster = false;
    odb::Rect sinksBbox = odb::Rect();
    sinksBbox.mergeInit();
    for (int branch : cluster) {
      odb::Point branchPt
          = {clkSteiner.branch[branch].x, clkSteiner.branch[branch].y};
      for (auto sink : point2pin[branchPt]) {
        if (!sink->isInputSignal()) {
          continue;
        }
        validCluster = true;
        int sinkX, sinkY;
        sink->getAvgXY(&sinkX, &sinkY);
        sinksBbox.merge({sinkX, sinkY});
      }
    }
    if (validCluster) {
      validClusters += 1;
      int64_t dist2Driver
          = odb::Point::manhattanDistance(sinksBbox.center(), driverPt);
      if (dist2Driver < minDist2Driver) {
        driverClusterID = n;
        minDist2Driver = dist2Driver;
      }
    }
  }

  debugPrint(logger_,
             CTS,
             "clock gate cloning",
             1,
             "Found {} clusters",
             validClusters);

  // Insert original ICG to its closest cluster, create clones to drive the
  // other clusters
  int nClones = 0;
  // hierarchy fix, make the clone net in the right scope
  sta::Pin* driver = nullptr;
  odb::dbModule* module
      = network_->getNetDriverParentModule(network_->dbToSta(icgNet), driver);
  if (module == nullptr) {
    // if none put in top level
    module = block_->getTopModule();
  }
  sta::Instance* scope
      = (module == nullptr || (module == block_->getTopModule()))
            ? network_->topInstance()
            : (sta::Instance*) (module->getModInst());

  for (int n = 0; n < clusters.size(); n++) {
    const std::vector<int>& cluster = clusters[n];
    if (cluster.empty()) {
      continue;
    }
    odb::dbInst* clone = nullptr;
    odb::dbNet* cloneNet = nullptr;
    bool disconectNets = true;
    odb::Rect sinksBbox = odb::Rect();
    sinksBbox.mergeInit();
    if (driverClusterID == n) {
      cloneNet = icgNet;
      clone = icgTerm->getInst();
      disconectNets = false;
      debugPrint(logger_,
                 CTS,
                 "clock gate cloning",
                 1,
                 "Original cell {}",
                 clone->getName());
      debugPrint(logger_,
                 CTS,
                 "clock gate cloning",
                 2,
                 " Original net {}",
                 cloneNet->getName());
    } else {
      // Create the ICG clone
      // Create a new input net
      std::string newNetName
          = "clonenet_" + std::to_string(++nClones) + "_" + icgNet->getName();

      cloneNet = network_->staToDb(network_->makeNet(
          newNetName.c_str(), scope, odb::dbNameUniquifyType::IF_NEEDED));

      cloneNet->setSigType(odb::dbSigType::CLOCK);

      staClockNets_.insert(cloneNet);

      // Create a new clone instance
      std::string newBufName
          = "clone_" + std::to_string(nClones) + "_" + icgName;
      odb::dbMaster* master = icgTerm->getInst()->getMaster();

      // fix: make buffer in same hierarchical module as driver
      clone = odb::dbInst::create(
          block_, master, newBufName.c_str(), false, module);

      clone->setSourceType(odb::dbSourceType::TIMING);

      debugPrint(logger_,
                 CTS,
                 "clock gate cloning",
                 1,
                 "Creating clone {} from {}",
                 newBufName,
                 icgName);
      debugPrint(logger_,
                 CTS,
                 "clock gate cloning",
                 2,
                 " New clone net {}",
                 cloneNet->getName());

      // Connect clone pins to same input nets as parent and new output net
      for (odb::dbITerm* iterm : clone->getITerms()) {
        if (iterm->isInputSignal()) {
          odb::dbITerm* parentITerm = icgTerm->getInst()->findITerm(
              iterm->getMTerm()->getName().c_str());
          odb::dbNet* parentNet = parentITerm->getNet();
          odb::dbModNet* parentModNet
              = network_->hierNet(network_->dbToSta(parentITerm));
          if (parentNet) {
            iterm->connect(parentNet);
            if (parentModNet) {
              iterm->connect(parentModNet);
            }
          }
        } else if (iterm->isOutputSignal()) {
          iterm->connect(cloneNet);
        }
      }
    }

    // Compute cluster center
    for (int branch : cluster) {
      odb::Point branchPt
          = {clkSteiner.branch[branch].x, clkSteiner.branch[branch].y};
      for (auto sink : point2pin[branchPt]) {
        if (!sink->isInputSignal()) {
          continue;
        }
        debugPrint(logger_,
                   CTS,
                   "clock gate cloning",
                   2,
                   "  Connects sink {}",
                   sink->getName());
        int sinkX, sinkY;
        sink->getAvgXY(&sinkX, &sinkY);
        sinksBbox.merge({sinkX, sinkY});
        if (disconectNets) {
          // Connect sinks to new clone instance
          sink->disconnect();
          sink->connect(cloneNet);
          sta::Pin* sinkPin = network_->dbToSta(sink);
          sta::Instance* sinkParentInst
              = network_->getOwningInstanceParent(sinkPin);
          if (sinkParentInst != scope) {
            network_->hierarchicalConnect(
                clone->getFirstOutput(), sink, cloneNet->getName().c_str());
          }
        }
      }
    }

    // Move ICG (clone or original) to the center of its sinks
    clone->setLocation(sinksBbox.xCenter(), sinksBbox.yCenter());
    clone->setPlacementStatus(odb::dbPlacementStatus::PLACED);
  }
  debugPrint(
      logger_, CTS, "clock gate cloning", 1, "Created {} clones", nClones);
}

void TritonCTS::populateTritonCTS()
{
  clearNumClocks();

  // Use dbSta to find all clock nets in the design.
  std::vector<std::pair<std::set<odb::dbNet*>, std::string>> clockNetsInfo;

  // Checks the user input in case there are other nets that need to be added to
  // the set.
  std::vector<odb::dbNet*> inputClkNets = options_->getClockNetsObjs();

  std::set<odb::dbNet*> allClkNets;
  if (!inputClkNets.empty()) {
    std::set<odb::dbNet*> clockNets;
    for (odb::dbNet* net : inputClkNets) {
      // Since a set is unique, only the nets not found by dbSta are added.
      clockNets.insert(net);
    }
    allClkNets.insert(clockNets.begin(), clockNets.end());
    clockNetsInfo.emplace_back(clockNets, "");
  } else {
    staClockNets_ = openSta_->findClkNets();
    sta::Sdc* sdc = openSta_->sdc();
    for (auto clk : *sdc->clocks()) {
      std::string clkName = clk->name();
      std::set<odb::dbNet*> clkNets;
      findClockRoots(clk, clkNets);
      for (auto net : clkNets) {
        if (allClkNets.find(net) != allClkNets.end()) {
          logger_->error(
              CTS, 114, "Clock {} overlaps a previous clock.", clkName);
        }
      }
      clockNetsInfo.emplace_back(clkNets, clkName);
      allClkNets.insert(clkNets.begin(), clkNets.end());
    }
  }
  // Iterate over all the nets found by the user-input and dbSta
  for (const auto& clockInfo : clockNetsInfo) {
    std::set<odb::dbNet*> clockNets = clockInfo.first;
    std::string clkName = clockInfo.second;
    for (odb::dbNet* net : clockNets) {
      if (net != nullptr) {
        cloneClockGaters(net);
        if (clkName.empty()) {
          logger_->info(CTS, 95, "Net \"{}\" found.", net->getName());
        } else {
          logger_->info(CTS,
                        7,
                        "Net \"{}\" found for clock \"{}\".",
                        net->getName(),
                        clkName);
        }
        // Initializes the net in TritonCTS. If the number of sinks is less than
        // 2, the net is discarded.
        if (visitedClockNets_.find(net) == visitedClockNets_.end()) {
          initOneClockTree(net, net, clkName, nullptr);
        }
      } else {
        logger_->warn(
            CTS,
            40,
            "Net was not found in the design for {}, please check. Skipping...",
            clkName);
      }
    }
  }

  if (getNumClocks() == 0) {
    logger_->warn(CTS, 83, "No clock nets have been found.");
  }

  logger_->info(CTS, 8, "TritonCTS found {} clock nets.", getNumClocks());
  options_->setNumClockRoots(getNumClocks());
}

TreeBuilder* TritonCTS::initClock(odb::dbNet* firstNet,
                                  odb::dbNet* clkInputNet,
                                  const std::string& sdcClock,
                                  TreeBuilder* parentBuilder)
{
  std::string driver;
  odb::dbITerm* iterm = firstNet->getFirstOutput();
  int xPin, yPin;
  if (iterm == nullptr) {
    odb::dbBTerm* bterm = firstNet->get1stBTerm();  // Clock pin
    if (bterm == nullptr) {
      logger_->info(
          CTS,
          122,
          "Clock net \"{}\" is skipped for CTS because it is not "
          "connected to any output instance pin or input block terminal.",
          firstNet->getName());
      return nullptr;
    }
    driver = bterm->getConstName();
    bterm->getFirstPinLocation(xPin, yPin);
  } else {
    odb::dbInst* inst = iterm->getInst();
    odb::dbMTerm* mterm = iterm->getMTerm();
    driver = std::string(inst->getConstName()) + "/"
             + std::string(mterm->getConstName());
    int xTmp, yTmp;
    computeITermPosition(iterm, xTmp, yTmp);
    xPin = xTmp;
    yPin = yTmp;
  }

  // Initialize clock net
  Clock clockNet(firstNet->getConstName(), driver, sdcClock, xPin, yPin);
  clockNet.setDriverPin(iterm);

  // Build a set of all the clock buffers' masters
  std::unordered_set<odb::dbMaster*> buffer_masters;
  for (const std::string& name : options_->getBufferList()) {
    auto master = db_->findMaster(name.c_str());
    if (master) {
      buffer_masters.insert(master);
    }
  }

  // Add the root buffer
  {
    const std::string& name = options_->getRootBuffer();
    auto master = db_->findMaster(name.c_str());
    if (master) {
      buffer_masters.insert(master);
    }
  }

  // Build a clock tree to drive macro cells with insertion delays
  // separated from registers or leaves without insertion delays
  TreeBuilder* builder = initClockTreeForMacrosAndRegs(
      firstNet, clkInputNet, buffer_masters, clockNet, parentBuilder);
  return builder;
}

// Build a separate clock tree to pull macro cells with insertion delays
// ahead of cells without insertion delays.  If sinks consist of
// both macros and FFs, clock tree for macros is built first.  A new net and a
// new buffer are created to drive cells without insertion delays.   New
// buffer will be sized later based on macro cell insertion delays.
//
//                |----|>----[] cells with insertion delays
//       firstNet |
//            |   |----|>----[]
//            v   |
//   [root]-------|                  |---|>----[] cells without insertion
//                |
//                |----|>------------|
//                      ^        ^   |
//                      |        |   |
//               new buffer secondNet|---|>----[]
//
TreeBuilder* TritonCTS::initClockTreeForMacrosAndRegs(
    odb::dbNet*& firstNet,
    odb::dbNet* clkInputNet,
    const std::unordered_set<odb::dbMaster*>& buffer_masters,
    Clock& clockNet,
    TreeBuilder* parentBuilder)
{
  // Separate sinks into two buckets: one with insertion delays and another
  // without
  std::vector<std::pair<odb::dbInst*, odb::dbMTerm*>> macroSinks;
  std::vector<std::pair<odb::dbInst*, odb::dbMTerm*>> registerSinks;
  if (!separateMacroRegSinks(
          firstNet, clockNet, buffer_masters, registerSinks, macroSinks)) {
    return nullptr;
  }

  // CTS_FORCE_MACRO_REG_SPLIT: force macro/register split even when
  // -no_insertion_delay is active. Without this, -no_insertion_delay
  // collapses everything into one unified tree (diverges from Pin3D).
  bool forceMacroRegSplit = false;
  if (const char* e = std::getenv("CTS_FORCE_MACRO_REG_SPLIT")) {
    forceMacroRegSplit = (std::atoi(e) != 0);
    if (forceMacroRegSplit && !macroSinks.empty() && !registerSinks.empty()) {
      logger_->info(CTS, 699,
          "CTS_FORCE_MACRO_REG_SPLIT: forcing split "
          "(macros={}, registers={})",
          macroSinks.size(), registerSinks.size());
    }
  }

  if ((!options_->insertionDelayEnabled() && !forceMacroRegSplit)
      || macroSinks.empty() || registerSinks.empty()) {
    // There is no need for separate clock trees
    for (odb::dbITerm* iterm : firstNet->getITerms()) {
      odb::dbInst* inst = iterm->getInst();
      if (iterm->isInputSignal() && inst->isPlaced()) {
        odb::dbMTerm* mterm = iterm->getMTerm();
        std::string name = std::string(inst->getConstName()) + "/"
                           + std::string(mterm->getConstName());
        int x, y;
        computeITermPosition(iterm, x, y);
        float insDelay = computeInsertionDelay(name, inst, mterm);
        clockNet.addSink(name, x, y, iterm, getInputPinCap(iterm), insDelay);
      }
    }
    if (clockNet.getNumSinks() < 2) {
      logger_->warn(CTS,
                    41,
                    "Net \"{}\" has {} sinks. Skipping...",
                    clockNet.getName(),
                    clockNet.getNumSinks());
      return nullptr;
    }
    logger_->info(CTS,
                  10,
                  " Clock net \"{}\" has {} sinks.",
                  firstNet->getConstName(),
                  clockNet.getNumSinks());
    int totalSinks = options_->getNumSinks() + clockNet.getNumSinks();
    options_->setNumSinks(totalSinks);
    incrementNumClocks();
    clockNet.setNetObj(firstNet);
    return addBuilder(
        options_, clockNet, clkInputNet, parentBuilder, logger_, db_);
  }

  // add macro sinks to existing firstNet
  TreeBuilder* firstBuilder = addClockSinks(
      clockNet, clkInputNet, firstNet, macroSinks, parentBuilder, "macros");
  if (firstBuilder) {
    firstBuilder->setTreeType(TreeType::MacroTree);
  }

  // create a new net 'secondNet' to drive register sinks
  odb::dbNet* secondNet;
  std::string topBufferName;
  Clock clockNet2 = forkRegisterClockNetwork(
      clockNet, registerSinks, firstNet, secondNet, topBufferName);

  // add register sinks to secondNet
  TreeBuilder* secondBuilder
      = addClockSinks(clockNet2,
                      clkInputNet,
                      secondNet,
                      registerSinks,
                      firstBuilder ? firstBuilder : parentBuilder,
                      "registers");
  if (secondBuilder) {
    secondBuilder->setTreeType(TreeType::RegisterTree);
    secondBuilder->setTopBufferName(std::move(topBufferName));
    secondBuilder->setDrivingNet(firstNet);
  }

  return firstBuilder;
}

// Separate sinks into registers (no insertion delay) and macros (insertion
// delay)
bool TritonCTS::separateMacroRegSinks(
    odb::dbNet*& net,
    Clock& clockNet,
    const std::unordered_set<odb::dbMaster*>& buffer_masters,
    std::vector<std::pair<odb::dbInst*, odb::dbMTerm*>>& registerSinks,
    std::vector<std::pair<odb::dbInst*, odb::dbMTerm*>>& macroSinks)
{
  for (odb::dbITerm* iterm : net->getITerms()) {
    odb::dbInst* inst = iterm->getInst();

    if (buffer_masters.find(inst->getMaster()) != buffer_masters.end()
        && inst->getSourceType() == odb::dbSourceType::TIMING) {
      logger_->warn(CTS,
                    105,
                    "Net \"{}\" already has clock buffer {}. Skipping...",
                    clockNet.getName(),
                    inst->getName());
      return false;
    }

    if (iterm->isInputSignal() && inst->isPlaced()) {
      // Cells with insertion delay, macros, clock gaters and inverters that
      // drive macros are put in the macro sinks.
      odb::dbMTerm* mterm = iterm->getMTerm();

      bool nonSinkMacro = !isSink(iterm);
      sta::Cell* masterCell = network_->dbToSta(mterm->getMaster());
      sta::LibertyCell* libertyCell = network_->libertyCell(masterCell);
      if (libertyCell && libertyCell->isInverter()) {
        odb::dbITerm* invertedTerm
            = inst->getFirstOutput()->getNet()->get1stSignalInput(false);
        nonSinkMacro &= isMacroBlockInst(invertedTerm->getInst());
      }

      // Macro bucket: real block macros (isBlock) on ANY tier +
      // bottom-tier name-matched macros only.
      // History: CODEX2 put all 28 macros → LatencyBalancer 226 delay bufs.
      // CODEX3 bottom-tier-only → 7 macros, but upper SRAMs (fakeram) in
      // register tree → 900ps latency at die corner → hold -320ps.
      // Fix: isBlock() SRAMs go to macro tree regardless of tier (shorter
      // tree, ~300ps). Name-matched-only macros still filtered by tier.
      // Safe because -no_insertion_delay disables LatencyBalancer.
      bool isMacro = hasInsertionDelay(inst, mterm) || nonSinkMacro
                     || isMacroBlockInst(inst);
      if (isMacro && cts3dDb_ != nullptr) {
        int instTier = cts3dDb_->getInstTier(inst);
        // Upper-tier: keep real block macros (inst->isBlock() or
        // master->isBlock()), filter out name-only matches.
        if (instTier > 0 && !inst->isBlock()
            && !(inst->getMaster() && inst->getMaster()->isBlock())) {
          isMacro = false;
        }
      }
      if (isMacro) {
        macroSinks.emplace_back(inst, mterm);
      } else {
        registerSinks.emplace_back(inst, mterm);
      }
    }
  }

  return true;
}

TreeBuilder* TritonCTS::addClockSinks(
    Clock& clockNet,
    odb::dbNet* topInputNet,
    odb::dbNet* physicalNet,
    const std::vector<std::pair<odb::dbInst*, odb::dbMTerm*>>& sinks,
    TreeBuilder* parentBuilder,
    const std::string& macrosOrRegs)
{
  for (auto elem : sinks) {
    odb::dbInst* inst = elem.first;
    odb::dbMTerm* mterm = elem.second;
    std::string name = std::string(inst->getConstName()) + "/"
                       + std::string(mterm->getConstName());
    int x, y;
    odb::dbITerm* iterm = inst->getITerm(mterm);
    computeITermPosition(iterm, x, y);
    float insDelay = computeInsertionDelay(name, inst, mterm);
    clockNet.addSink(name, x, y, iterm, getInputPinCap(iterm), insDelay);
  }
  logger_->info(CTS,
                11,
                " Clock net \"{}\" for {} has {} sinks.",
                physicalNet->getConstName(),
                macrosOrRegs,
                clockNet.getNumSinks());
  int totalSinks = options_->getNumSinks() + clockNet.getNumSinks();
  options_->setNumSinks(totalSinks);
  incrementNumClocks();
  clockNet.setNetObj(physicalNet);
  return addBuilder(
      options_, clockNet, topInputNet, parentBuilder, logger_, db_);
}

Clock TritonCTS::forkRegisterClockNetwork(
    Clock& clockNet,
    const std::vector<std::pair<odb::dbInst*, odb::dbMTerm*>>& registerSinks,
    odb::dbNet*& firstNet,
    odb::dbNet*& secondNet,
    std::string& topBufferName)
{
  // create a new clock net to drive register sinks
  std::string newClockName = clockNet.getName() + "_regs";
  secondNet = odb::dbNet::create(block_, newClockName.c_str());
  secondNet->setSigType(odb::dbSigType::CLOCK);

  sta::Pin* first_pin_driver = nullptr;
  odb::dbModule* first_net_module = network_->getNetDriverParentModule(
      network_->dbToSta(firstNet), first_pin_driver);
  (void) first_pin_driver;
  sta::Pin* second_pin_driver = nullptr;
  odb::dbModule* second_net_module = network_->getNetDriverParentModule(
      network_->dbToSta(secondNet), second_pin_driver);
  (void) second_pin_driver;
  odb::dbModule* target_module = nullptr;
  if ((first_net_module != nullptr)
      && (first_net_module == second_net_module)) {
    target_module = first_net_module;
  }

  // move register sinks from previous clock net to new clock net
  for (auto elem : registerSinks) {
    odb::dbInst* inst = elem.first;
    odb::dbMTerm* mterm = elem.second;
    odb::dbITerm* iterm = inst->getITerm(mterm);
    iterm->disconnect();
    iterm->connect(secondNet);
  }

  // create a new clock buffer
  odb::dbMaster* master = db_->findMaster(options_->getRootBuffer().c_str());
  topBufferName = "clkbuf_regs_" + std::to_string(regTreeRootBufIndex_++) + "_"
                  + clockNet.getSdcName();
  odb::dbInst* clockBuf = odb::dbInst::create(
      block_, master, topBufferName.c_str(), false, target_module);

  // place new clock buffer near center of mass for registers
  odb::Rect bbox = secondNet->getTermBBox();
  clockBuf->setSourceType(odb::dbSourceType::TIMING);
  clockBuf->setLocation(bbox.xCenter(), bbox.yCenter());
  clockBuf->setPlacementStatus(odb::dbPlacementStatus::PLACED);

  // connect root buffer to clock net
  odb::dbITerm* inputTerm = getFirstInput(clockBuf);
  odb::dbITerm* outputTerm = clockBuf->getFirstOutput();
  inputTerm->connect(firstNet);
  outputTerm->connect(secondNet);

  // initialize new clock net
  std::string driver = std::string(clockBuf->getConstName()) + "/"
                       + std::string(outputTerm->getMTerm()->getConstName());
  int xPin, yPin;
  computeITermPosition(outputTerm, xPin, yPin);
  Clock clockNet2(
      secondNet->getConstName(), driver, clockNet.getSdcName(), xPin, yPin);
  clockNet2.setDriverPin(outputTerm);

  return clockNet2;
}

void TritonCTS::computeITermPosition(odb::dbITerm* term, int& x, int& y) const
{
  odb::dbITermShapeItr itr;

  odb::dbShape shape;
  x = 0;
  y = 0;
  unsigned numShapes = 0;
  for (itr.begin(term); itr.next(shape);) {
    if (!shape.isVia()) {
      x += shape.xMin() + (shape.xMax() - shape.xMin()) / 2;
      y += shape.yMin() + (shape.yMax() - shape.yMin()) / 2;
      ++numShapes;
    }
  }
  if (numShapes > 0) {
    x /= numShapes;
    y /= numShapes;
  }
};

void TritonCTS::destroyClockModNet(sta::Pin* pin_driver)
{
  if (pin_driver == nullptr || network_->hasHierarchy() == false) {
    return;
  }

  odb::dbModNet* mod_net = network_->hierNet(pin_driver);
  if (mod_net) {
    odb::dbModNet::destroy(mod_net);
  }
}

void TritonCTS::writeClockNetsToDb(TreeBuilder* builder,
                                   std::set<odb::dbNet*>& clkLeafNets)
{
  Clock& clockNet = builder->getClock();
  odb::dbNet* topClockNet = clockNet.getNetObj();

  logger_->info(CTS, 880,
      "writeClockNetsToDb: clock='{}', topNet='{}', treeType={}, "
      "topBufferName='{}'",
      clockNet.getName(),
      topClockNet ? topClockNet->getConstName() : "NULL",
      static_cast<int>(builder->getTreeType()),
      builder->getTopBufferName());

  // gets the module for the driver for the net
  sta::Pin* pin_driver = nullptr;
  odb::dbModule* top_module = network_->getNetDriverParentModule(
      network_->dbToSta(topClockNet), pin_driver);
  (void) pin_driver;

  disconnectAllSinksFromNet(topClockNet);
  destroyClockModNet(pin_driver);

  // re-connect top buffer that separates macros from registers
  if (builder->getTreeType() == TreeType::RegisterTree) {
    odb::dbInst* topRegBuffer
        = block_->findInst(builder->getTopBufferName().c_str());
    if (topRegBuffer) {
      odb::dbITerm* topRegBufferInputPin = getFirstInput(topRegBuffer);
      topRegBufferInputPin->connect(builder->getDrivingNet());
      logger_->info(CTS, 881,
          "  RegisterTree: connected '{}' input to drivingNet '{}'",
          builder->getTopBufferName(),
          builder->getDrivingNet() ? builder->getDrivingNet()->getConstName() : "NULL");
    } else {
      logger_->warn(CTS, 882,
          "  RegisterTree: topRegBuffer '{}' NOT FOUND in ODB!",
          builder->getTopBufferName());
    }
  }

  createClockBuffers(clockNet, top_module);

  // Connect ALL root buffers to top clock net.
  {
    const std::string clkName = clockNet.getName();
    std::vector<std::string> rootBufCandidates;
    rootBufCandidates.push_back("clkbuf_0_" + clkName);
    rootBufCandidates.push_back("t1_clkbuf_0_" + clkName);
    bool anyConnected = false;
    for (const auto& rootName : rootBufCandidates) {
      odb::dbInst* rootInst = block_->findInst(rootName.c_str());
      if (rootInst) {
        odb::dbITerm* rootInputPin = getFirstInput(rootInst);
        if (rootInputPin) {
          rootInputPin->connect(topClockNet);
          anyConnected = true;
          logger_->info(CTS, 883,
              "  Root '{}' input connected to topClockNet '{}'",
              rootName, topClockNet->getConstName());
        }
      } else {
        logger_->info(CTS, 884,
            "  Root candidate '{}' not found in ODB", rootName);
      }
    }
    if (!anyConnected) {
      logger_->error(CTS, 498,
          "No root clock buffer found for clock '{}'", clkName);
    }
  }
  topClockNet->setSigType(odb::dbSigType::CLOCK);

  std::map<int, int> fanoutcount;

  // create subNets
  numClkNets_ = 0;
  numFixedNets_ = 0;
  ClockSubNet* rootSubNet = nullptr;
  std::set<ClockInst*> removedSinks;
  clockNet.forEachSubNet([&](ClockSubNet& subNet) {
    bool outputPinFound = true;
    bool inputPinFound = true;
    bool leafLevelNet = subNet.isLeafLevel();
    // Match root subnet with or without tier prefix
    const std::string subNetName = subNet.getName();
    const std::string clkName2 = clockNet.getName();
    if (subNetName == ("clknet_0_" + clkName2)
        || subNetName == ("t1_clknet_0_" + clkName2)) {
      rootSubNet = &subNet;
    }

    // Fix A: Early skip for orphan subnets (numSinks=0 and not root).
    // CTS-0080 float precision skip can leave leaf clusters with 0 sinks.
    // Their driver buffer exists in ODB but has no connected sinks.
    // Proceeding would create a disconnected ODB net and risk SIGSEGV
    // when accessing the driver's input pin (no upstream net yet).
    // Also mark the orphan driver in removedSinks so branchBufferCount()
    // skips it during recursive tree traversal (otherwise it reaches the
    // orphan via parent subnet → outITerm->getNet() = nullptr → SIGSEGV).
    if (subNet.getNumSinks() == 0 && &subNet != rootSubNet) {
      logger_->warn(CTS, 715,
          "  Orphan subnet '{}' (0 sinks, driver='{}') — skipping ODB net creation",
          subNet.getName(),
          subNet.getDriver() ? subNet.getDriver()->getName() : "NULL");
      if (subNet.getDriver()) {
        removedSinks.insert(subNet.getDriver());
      }
      return;  // lambda early return — skip this subnet entirely
    }

    odb::dbNet* clkSubNet
        = odb::dbNet::create(block_, subNet.getName().c_str());
    subNet.setNetObj(clkSubNet);

    ++numClkNets_;
    clkSubNet->setSigType(odb::dbSigType::CLOCK);

    odb::dbInst* driver = subNet.getDriver()->getDbInst();
    odb::dbITerm* driverInputPin = getFirstInput(driver);
    odb::dbNet* inputNet = driverInputPin->getNet();
    odb::dbITerm* outputPin = driver->getFirstOutput();
    if (outputPin == nullptr) {
      outputPinFound = false;
    }
    if (outputPinFound) {
      outputPin->connect(clkSubNet);
    }

    if (subNet.getNumSinks() == 0) {
      inputPinFound = false;
    }

    subNet.forEachSink([&](ClockInst* inst) {
      // Skip orphan drivers — their subnet was skipped (CTS-0715).
      // Without this, orphan buffer's input gets connected to parent net,
      // allowing branchBufferCount to reach it → dangling output net → SIGSEGV.
      if (removedSinks.find(inst) != removedSinks.end()) {
        return;
      }
      odb::dbITerm* inputPin = nullptr;
      if (inst->isClockBuffer()) {
        odb::dbInst* sink = inst->getDbInst();
        inputPin = getFirstInput(sink);
      } else {
        inputPin = inst->getDbInputPin();
      }
      if (inputPin == nullptr) {
        inputPinFound = false;
      } else {
        if (!inputPin->getInst()->isPlaced()) {
          inputPinFound = false;
        }
      }

      if (inputPinFound) {
        inputPin->connect(clkSubNet);
        // get module for input pin
        // resolve connection in hierarchy
        if (network_->hasHierarchy()) {
          network_->hierarchicalConnect(
              outputPin, inputPin, clkSubNet->getName().c_str());
        }
      }
    });

    if (leafLevelNet) {
      // Report fanout values only for sink nets
      if (fanoutcount.find(subNet.getNumSinks()) == fanoutcount.end()) {
        fanoutcount[subNet.getNumSinks()] = 0;
      }
      fanoutcount[subNet.getNumSinks()] = fanoutcount[subNet.getNumSinks()] + 1;
      clkLeafNets.insert(clkSubNet);
    }

    if (!inputPinFound || !outputPinFound) {
      // Net not fully connected. Removing it.
      logger_->warn(CTS, 885,
          "  Disconnected subnet '{}' (driver='{}', "
          "inputPinFound={}, outputPinFound={}, numSinks={}) — keeping as-is",
          subNet.getName(),
          driver ? driver->getConstName() : "NULL",
          inputPinFound, outputPinFound, subNet.getNumSinks());
      // V58: Do NOT destroy disconnected subnets — destroy cascade can cause
      // SIGSEGV in checkUpstreamConnections when upstream nets/insts are
      // already freed. Orphan nets are harmless (no timing impact, DRT skips them).
    }
  });

  if (!rootSubNet) {
    logger_->error(
        CTS, 85, "Could not find the root of {}", clockNet.getName());
  }

  int minPath = std::numeric_limits<int>::max();
  int maxPath = std::numeric_limits<int>::min();
  rootSubNet->forEachSink([&](ClockInst* inst) {
    // skip removed sinks
    if (removedSinks.find(inst) == removedSinks.end()) {
      if (inst->isClockBuffer()) {
        std::pair<int, int> resultsForBranch
            = branchBufferCount(inst, 1, clockNet);
        minPath = std::min(resultsForBranch.first, minPath);
        maxPath = std::max(resultsForBranch.second, maxPath);
      }
    } else {
      rootSubNet->removeSinks(removedSinks);
    }
  });

  logger_->info(
      CTS, 12, "    Minimum number of buffers in the clock path: {}.", minPath);
  logger_->info(
      CTS, 13, "    Maximum number of buffers in the clock path: {}.", maxPath);

  if (numFixedNets_ > 0) {
    logger_->info(
        CTS, 14, "    {} clock nets were removed/fixed.", numFixedNets_);
  }

  logger_->info(CTS, 15, "    Created {} clock nets.", numClkNets_);
  int totalNets = options_->getNumClockSubnets() + numClkNets_;
  options_->setNumClockSubnets(totalNets);

  std::string fanout;
  for (auto const& x : fanoutcount) {
    fanout += std::to_string(x.first) + ':' + std::to_string(x.second) + ", ";
  }

  logger_->info(CTS,
                16,
                "    Fanout distribution for the current clock = {}.",
                fanout.substr(0, fanout.size() - 2) + ".");
  logger_->info(
      CTS, 17, "    Max level of the clock tree: {}.", clockNet.getMaxLevel());
}

// Utility function to get all unique clock tree levels
std::vector<int> TritonCTS::getAllClockTreeLevels(Clock& clockNet)
{
  std::set<int> uniqueLevels;

  clockNet.forEachSubNet([&](ClockSubNet& subNet) {
    if (!subNet.isLeafLevel() && subNet.getTreeLevel() != -1) {
      uniqueLevels.insert(subNet.getTreeLevel());
    }
  });

  return std::vector<int>(uniqueLevels.begin(), uniqueLevels.end());
}

// Function to apply NDR to specific clock tree levels and return the number of
// NDR applied nets
int TritonCTS::applyNDRToClockLevels(Clock& clockNet,
                                     odb::dbTechNonDefaultRule* clockNDR,
                                     const std::vector<int>& targetLevels)
{
  int ndrAppliedNets = 0;

  debugPrint(
      logger_, CTS, "clustering", 1, "Applying NDR to clock tree levels: ");
  for (int level : targetLevels) {
    debugPrint(logger_, CTS, "clustering", 1, "{} ", level);
  }

  // Check if the main clock net (level 0) is in the level list
  if (std::ranges::find(targetLevels, 0) != targetLevels.end()) {
    odb::dbNet* clk_net = clockNet.getNetObj();
    clk_net->setNonDefaultRule(clockNDR);
    ndrAppliedNets++;
    // clang-format off
    debugPrint(logger_, CTS, "clustering", 1,
        "Applied NDR to: {} (level {})", clockNet.getName(), 0);
    // clang-format on
  }

  // Check clock sub nets list and apply NDR if level matches
  clockNet.forEachSubNet([&](ClockSubNet& subNet) {
    int level = subNet.getTreeLevel();
    if (std::ranges::find(targetLevels, level) != targetLevels.end()) {
      odb::dbNet* net = subNet.getNetObj();
      if (!subNet.isLeafLevel()) {
        net->setNonDefaultRule(clockNDR);
        ndrAppliedNets++;
        std::string net_name = net->getName();
        // clang-format off
        debugPrint(logger_, CTS, "clustering", 1,
            "Applied NDR to: {} (level {})", net_name, level);
        // clang-format on
      }
    }
  });

  return ndrAppliedNets;
}

// Alternative function to apply NDR to a range of clock tree levels
int TritonCTS::applyNDRToClockLevelRange(Clock& clockNet,
                                         odb::dbTechNonDefaultRule* clockNDR,
                                         const int minLevel,
                                         const int maxLevel)
{
  std::vector<int> targetLevels;
  for (int i = minLevel; i <= maxLevel; i++) {
    targetLevels.push_back(i);
  }

  return applyNDRToClockLevels(clockNet, clockNDR, targetLevels);
}

// Function to apply NDR to the first half of clock tree levels
int TritonCTS::applyNDRToFirstHalfLevels(Clock& clockNet,
                                         odb::dbTechNonDefaultRule* clockNDR)
{
  // Get all unique levels in the design
  const std::vector<int> allLevels = getAllClockTreeLevels(clockNet);

  // Calculate first half (rounding up if odd number of levels)
  const int halfCount = (allLevels.size() + 1) / 2;

  // Create vector with first half of levels
  std::vector<int> firstHalfLevels(allLevels.begin(),
                                   allLevels.begin() + halfCount);

  // clang-format off
  debugPrint(logger_, CTS, "clustering", 1, "Total clock tree levels found: {}"
        " Applying NDR to first {} levels", allLevels.size(), halfCount);
  // clang-format on

  // Apply NDR to the first half
  return applyNDRToClockLevels(clockNet, clockNDR, firstHalfLevels);
}

// Priority for minSpc rule is SPACINGTABLE TWOWIDTHS > SPACINGTABLE PRL >
// SPACING
int TritonCTS::getNetSpacing(odb::dbTechLayer* layer,
                             const int width1,
                             const int width2)
{
  int min_spc = 0;
  if (layer->hasTwoWidthsSpacingRules()) {
    min_spc = layer->findTwSpacing(width1, width2, 0);
  } else if (layer->hasV55SpacingRules()) {
    min_spc = layer->findV55Spacing(std::max(width1, width2), 0);
  } else if (!layer->getV54SpacingRules().empty()) {
    for (auto rule : layer->getV54SpacingRules()) {
      if (rule->hasRange()) {
        uint32_t rmin;
        uint32_t rmax;
        rule->getRange(rmin, rmax);
        if (width1 < rmin || width2 > rmax) {
          continue;
        }
      }
      min_spc = std::max<int>(min_spc, rule->getSpacing());
    }
  } else {
    min_spc = layer->getSpacing();
  }

  // Last resort, get pitch - minWidth
  if (min_spc == 0) {
    min_spc = layer->getPitch() - layer->getMinWidth();
  }

  return min_spc;
}

void TritonCTS::writeClockNDRsToDb(TreeBuilder* builder)
{
  char ruleName[64];
  int ruleIndex = 0;
  odb::dbTechNonDefaultRule* clockNDR;
  Clock& clockNet = builder->getClock();

  // create a new non-default rule in *block* not tech
  while (ruleIndex >= 0) {
    snprintf(ruleName, 64, "CTS_NDR_%d", ruleIndex++);
    clockNDR = odb::dbTechNonDefaultRule::create(block_, ruleName);
    if (clockNDR) {
      break;
    }
  }
  assert(clockNDR != nullptr);

  // define NDR for all routing layers
  odb::dbTech* tech = db_->getTech();
  for (int i = 1; i <= tech->getRoutingLayerCount(); i++) {
    odb::dbTechLayer* layer = tech->findRoutingLayer(i);
    odb::dbTechLayerRule* layerRule = clockNDR->getLayerRule(layer);
    if (!layerRule) {
      layerRule = odb::dbTechLayerRule::create(clockNDR, layer);
    }
    assert(layerRule != nullptr);

    int defaultWidth = layer->getWidth();
    int defaultSpace = getNetSpacing(layer, defaultWidth, defaultWidth);

    // If width or space is 0, something is not right
    if (defaultWidth == 0 || defaultSpace == 0) {
      logger_->warn(CTS,
                    208,
                    "Clock NDR settings for layer {}: defaultSpace: {} - "
                    "defaultWidth: {}",
                    layer->getName(),
                    defaultSpace,
                    defaultWidth);
    }

    // Set NDR settings
    int ndr_width = defaultWidth;
    layerRule->setWidth(ndr_width);
    int ndr_space = 2 * getNetSpacing(layer, ndr_width, ndr_width);
    layerRule->setSpacing(ndr_space);

    debugPrint(logger_,
               CTS,
               "clustering",
               1,
               "  NDR rule set to layer {} {} as space={} width={} vs. default "
               "space={} width={}",
               i,
               layer->getName(),
               layerRule->getSpacing(),
               layerRule->getWidth(),
               defaultSpace,
               defaultWidth);
  }

  int clkNets = 0;

  // Apply NDR following the selected strategy (root_only, half, full)
  switch (options_->getApplyNdr()) {
    case CtsOptions::NdrStrategy::ROOT_ONLY:
      clkNets = applyNDRToClockLevels(clockNet, clockNDR, {0});
      break;
    case CtsOptions::NdrStrategy::HALF:
      clkNets = applyNDRToFirstHalfLevels(clockNet, clockNDR);
      break;
    case CtsOptions::NdrStrategy::FULL:
      clkNets = applyNDRToClockLevels(
          clockNet, clockNDR, getAllClockTreeLevels(clockNet));
      break;
    case CtsOptions::NdrStrategy::NONE:
      // Should not be called
      break;
  }

  debugPrint(logger_,
             CTS,
             "clustering",
             1,
             "Non-default rule {} for double spacing has been applied to {} "
             "clock nets",
             ruleName,
             clkNets);
}

std::pair<int, int> TritonCTS::branchBufferCount(ClockInst* inst,
                                                 int bufCounter,
                                                 Clock& clockNet)
{
  odb::dbInst* sink = inst->getDbInst();
  // Null guard: orphan buffers (Fix A skip) have no ODB net on output.
  // Return current count as both min/max to avoid SIGSEGV.
  if (sink == nullptr) {
    return {bufCounter, bufCounter};
  }
  odb::dbITerm* outITerm = sink->getFirstOutput();
  if (outITerm == nullptr || outITerm->getNet() == nullptr) {
    return {bufCounter, bufCounter};
  }
  int minPath = std::numeric_limits<int>::max();
  int maxPath = std::numeric_limits<int>::min();
  for (odb::dbITerm* sinkITerms : outITerm->getNet()->getITerms()) {
    if (sinkITerms != outITerm) {
      ClockInst* clockInst
          = clockNet.findClockByName(sinkITerms->getInst()->getName());
      if (clockInst == nullptr) {
        int newResult = bufCounter + 1;
        maxPath = std::max(newResult, maxPath);
        minPath = std::min(newResult, minPath);
      } else {
        std::pair<int, int> newResults
            = branchBufferCount(clockInst, bufCounter + 1, clockNet);
        minPath = std::min(newResults.first, minPath);
        maxPath = std::max(newResults.second, maxPath);
      }
    }
  }
  std::pair<int, int> results(minPath, maxPath);
  return results;
}

void TritonCTS::disconnectAllSinksFromNet(odb::dbNet* net)
{
  odb::dbSet<odb::dbITerm> iterms = net->getITerms();
  for (odb::dbITerm* iterm : iterms) {
    if (iterm->getIoType() == odb::dbIoType::INPUT) {
      iterm->disconnect();
    }
  }
}

void TritonCTS::disconnectAllPinsFromNet(odb::dbNet* net)
{
  odb::dbSet<odb::dbITerm> iterms = net->getITerms();
  for (odb::dbITerm* iterm : iterms) {
    iterm->disconnect();
  }
}

void TritonCTS::checkUpstreamConnections(odb::dbNet* net)
{
  while (net->getITermCount() <= 1) {
    // Net is incomplete, only 1 pin.
    odb::dbITerm* firstITerm = net->get1stITerm();
    if (firstITerm == nullptr) {
      disconnectAllPinsFromNet(net);
      odb::dbNet::destroy(net);
      break;
    }
    odb::dbInst* bufferInst = firstITerm->getInst();
    odb::dbITerm* driverInputPin = getFirstInput(bufferInst);
    disconnectAllPinsFromNet(net);
    odb::dbNet::destroy(net);
    net = driverInputPin->getNet();
    ++numFixedNets_;
    --numClkNets_;
    odb::dbInst::destroy(bufferInst);
  }
}

void TritonCTS::createClockBuffers(Clock& clockNet, odb::dbModule* parent)
{
  // scan for duplicate names in clockBuffers_ deque
  {
    std::unordered_map<std::string, int> nameCounts;
    clockNet.forEachClockBuffer([&](const ClockInst& inst) {
      nameCounts[inst.getName()]++;
    });
    for (const auto& [name, count] : nameCounts) {
      if (count > 1) {
        logger_->warn(CTS, 497,
            "Diag: buffer '{}' appears {} times in clockBuffers_ deque",
            name, count);
      }
    }
    // Also check if any buffer name already exists as ODB instance
    clockNet.forEachClockBuffer([&](const ClockInst& inst) {
      odb::dbInst* existing = block_->findInst(inst.getName().c_str());
      if (existing) {
        logger_->warn(CTS, 496,
            "Diag: buffer '{}' already exists in ODB BEFORE createClockBuffers "
            "(master={})",
            inst.getName(), existing->getMaster()->getName());
      }
    });
  }

  unsigned numBuffers = 0;
  clockNet.forEachClockBuffer([&](ClockInst& inst) {
    odb::dbMaster* master = db_->findMaster(inst.getMaster().c_str());
    odb::dbInst* newInst = odb::dbInst::create(
        block_, master, inst.getName().c_str(), false, parent);
    // null check to catch duplicate buffer names (dbInst::create
    // returns nullptr when an instance with the same name already exists).
    if (!newInst) {
      logger_->error(CTS, 499,
                     "Failed to create clock buffer '{}' — duplicate name?",
                     inst.getName());
    }
    newInst->setSourceType(odb::dbSourceType::TIMING);
    inst.setInstObj(newInst);
    inst2clkbuf_[newInst] = &inst;
    inst.setInputPinObj(getFirstInput(newInst));
    newInst->setLocation(inst.getX(), inst.getY());
    newInst->setPlacementStatus(odb::dbPlacementStatus::PLACED);
    ++numBuffers;
  });
  logger_->info(CTS, 18, "    Created {} clock buffers.", numBuffers);
  int totalBuffers = options_->getNumBuffersInserted() + numBuffers;
  options_->setNumBuffersInserted(totalBuffers);
}

odb::dbITerm* TritonCTS::getFirstInput(odb::dbInst* inst) const
{
  odb::dbSet<odb::dbITerm> iterms = inst->getITerms();
  for (odb::dbITerm* iterm : iterms) {
    if (iterm->isInputSignal()) {
      return iterm;
    }
  }

  return nullptr;
}

odb::dbITerm* TritonCTS::getSingleOutput(odb::dbInst* inst,
                                         odb::dbITerm* input) const
{
  odb::dbSet<odb::dbITerm> iterms = inst->getITerms();
  odb::dbITerm* output = nullptr;
  for (odb::dbITerm* iterm : iterms) {
    if (iterm != input && iterm->isOutputSignal()) {
      odb::dbNet* net = iterm->getNet();
      if (net) {
        if (staClockNets_.find(net) != staClockNets_.end()) {
          output = iterm;
          break;
        }
      }
    }
  }
  return output;
}
bool TritonCTS::masterExists(const std::string& master) const
{
  return db_->findMaster(master.c_str());
};

void TritonCTS::findClockRoots(sta::Clock* clk,
                               std::set<odb::dbNet*>& clockNets)
{
  std::vector<odb::dbNet*> skipNets = options_->getSkipNets();
  for (const sta::Pin* pin : clk->leafPins()) {
    odb::dbITerm* instTerm;
    odb::dbBTerm* port;
    odb::dbModITerm* moditerm;
    network_->staToDb(pin, instTerm, port, moditerm);
    odb::dbNet* net = instTerm ? instTerm->getNet() : port->getNet();
    if (std::ranges::find(skipNets, net) != skipNets.end()) {
      logger_->warn(CTS,
                    42,
                    "Skipping root net {}, specified by the user...",
                    net->getName());
      continue;
    }
    clockNets.insert(net);
  }
}

float TritonCTS::getInputPinCap(odb::dbITerm* iterm)
{
  odb::dbInst* inst = iterm->getInst();
  sta::Cell* masterCell = network_->dbToSta(inst->getMaster());
  sta::LibertyCell* libertyCell = network_->libertyCell(masterCell);
  if (!libertyCell) {
    return 0.0;
  }

  sta::LibertyPort* inputPort
      = libertyCell->findLibertyPort(iterm->getMTerm()->getConstName());
  if (inputPort) {
    return inputPort->capacitance();
  }

  return 0.0;
}

bool TritonCTS::isSink(odb::dbITerm* iterm)
{
  odb::dbInst* inst = iterm->getInst();
  sta::Cell* masterCell = network_->dbToSta(inst->getMaster());
  sta::LibertyCell* libertyCell = network_->libertyCell(masterCell);
  if (!libertyCell) {
    return true;
  }

  if (isMacroBlockInst(inst)) {
    return true;
  }

  sta::LibertyPort* inputPort
      = libertyCell->findLibertyPort(iterm->getMTerm()->getConstName());
  if (inputPort) {
    return inputPort->isRegClk();
  }

  return false;
}

bool TritonCTS::hasInsertionDelay(odb::dbInst* inst, odb::dbMTerm* mterm)
{
  if (options_->insertionDelayEnabled()) {
    sta::LibertyCell* libCell = network_->libertyCell(network_->dbToSta(inst));
    if (libCell) {
      sta::LibertyPort* libPort
          = libCell->findLibertyPort(mterm->getConstName());
      if (libPort) {
        const float rise = libPort->clkTreeDelay(
            0.0, sta::RiseFall::rise(), sta::MinMax::max());
        const float fall = libPort->clkTreeDelay(
            0.0, sta::RiseFall::fall(), sta::MinMax::max());

        if (rise != 0 || fall != 0) {
          return true;
        }
      }
    }
  }
  return false;
}

double TritonCTS::computeInsertionDelay(const std::string& name,
                                        odb::dbInst* inst,
                                        odb::dbMTerm* mterm)
{
  double insDelayPerMicron = 0.0;

  if (!options_->insertionDelayEnabled()) {
    return insDelayPerMicron;
  }

  sta::LibertyCell* libCell = network_->libertyCell(network_->dbToSta(inst));
  if (libCell) {
    sta::LibertyPort* libPort = libCell->findLibertyPort(mterm->getConstName());
    if (libPort) {
      const float rise = libPort->clkTreeDelay(
          0.0, sta::RiseFall::rise(), sta::MinMax::max());
      const float fall = libPort->clkTreeDelay(
          0.0, sta::RiseFall::fall(), sta::MinMax::max());

      if (rise != 0 || fall != 0) {
        // use average of max rise and max fall
        // TODO: do we need to look at min insertion delays?
        double delayPerSec = (rise + fall) / 2.0;
        // convert delay to length because HTree uses lengths
        sta::Corner* corner = openSta_->cmdCorner();
        double capPerMicron
            = estimate_parasitics_->wireSignalCapacitance(corner) * 1e-6;
        double resPerMicron
            = estimate_parasitics_->wireSignalResistance(corner) * 1e-6;
        if (sta::fuzzyEqual(capPerMicron, 1e-18)
            || sta::fuzzyEqual(resPerMicron, 1e-18)) {
          logger_->warn(CTS,
                        203,
                        "Insertion delay cannot be used because unit "
                        "capacitance or unit resistance is zero.  Check "
                        "layer RC settings.");
          return 0.0;
        }
        insDelayPerMicron = delayPerSec / (capPerMicron * resPerMicron);
        // clang-format off
          debugPrint(logger_, CTS, "clustering", 1, "sink {} has ins "
                     "delay={:.2e} and micron leng={:0.1f} dbUnits/um={}",
                     name, delayPerSec, insDelayPerMicron,
                     block_->getDbUnitsPerMicron());
          debugPrint(logger_, CTS, "clustering", 1, "capPerMicron={:.2e} "
                     "resPerMicron={:.2e}", capPerMicron, resPerMicron);
        // clang-format on
      }
    }
  }
  return insDelayPerMicron;
}

static float getInputCap(const sta::LibertyCell* cell)
{
  sta::LibertyPort *in, *out;
  cell->bufferPorts(in, out);
  if (in != nullptr) {
    return in->capacitance();
  }
  return 0.0;
}

static sta::LibertyCell* findBestDummyCell(
    const std::vector<sta::LibertyCell*>& dummyCandidates,
    float deltaCap)
{
  float minDiff = std::numeric_limits<float>::max();
  sta::LibertyCell* bestCell = nullptr;
  for (sta::LibertyCell* cell : dummyCandidates) {
    float diff = std::abs(getInputCap(cell) - deltaCap);
    if (diff < minDiff) {
      minDiff = diff;
      bestCell = cell;
    }
  }
  return bestCell;
}

int TritonCTS::writeDummyLoadsToDb(Clock& clockNet,
                                   std::unordered_set<odb::dbInst*>& dummies)
{
  // Traverse clock tree and compute ideal output caps for clock
  // buffers in the same level
  if (!computeIdealOutputCaps(clockNet)) {
    // No cap adjustment is needed
    return 0;
  }

  // Find suitable candidate cells for dummy loads
  std::vector<sta::LibertyCell*> dummyCandidates;
  findCandidateDummyCells(dummyCandidates);

  int nDummies = 0;
  clockNet.forEachSubNet([&](ClockSubNet& subNet) {
    subNet.forEachSink([&](ClockInst* inst) {
      if (inst->isClockBuffer()
          && !sta::fuzzyEqual(inst->getOutputCap(),
                              inst->getIdealOutputCap())) {
        odb::dbInst* dummyInst
            = insertDummyCell(clockNet, inst, dummyCandidates);
        if (dummyInst != nullptr) {
          dummies.insert(dummyInst);
          nDummies++;
        }
      }
    });
  });

  if (logger_->debugCheck(utl::CTS, "dummy load", 1)) {
    printClockNetwork(clockNet);
  }
  return nDummies;
}

// Return true if any clock buffers need cap adjustment; false otherwise
bool TritonCTS::computeIdealOutputCaps(Clock& clockNet)
{
  bool needAdjust = false;

  // pass 1: compute actual output caps seen by each clock instance
  clockNet.forEachSubNet([&](ClockSubNet& subNet) {
    // build driver -> subNet map
    ClockInst* driver = subNet.getDriver();
    driver2subnet_[driver] = &subNet;
    float sinkCapTotal = 0.0;
    subNet.forEachSink([&](ClockInst* inst) {
      odb::dbITerm* inputPin = inst->isClockBuffer()
                                   ? getFirstInput(inst->getDbInst())
                                   : inst->getDbInputPin();
      float cap = getInputPinCap(inputPin);
      // TODO: include wire caps?
      sinkCapTotal += cap;
    });
    driver->setOutputCap(sinkCapTotal);
  });

  // pass 2: compute ideal output caps for perfectly balanced tree
  clockNet.forEachSubNet([&](const ClockSubNet& subNet) {
    ClockInst* driver = subNet.getDriver();
    float maxCap = std::numeric_limits<float>::min();
    subNet.forEachSink([&](ClockInst* inst) {
      if (inst->isClockBuffer() && inst->getOutputCap() > maxCap) {
        maxCap = inst->getOutputCap();
      }
    });
    subNet.forEachSink([&](ClockInst* inst) {
      if (inst->isClockBuffer()) {
        inst->setIdealOutputCap(maxCap);
        float cap = inst->getOutputCap();
        if (!sta::fuzzyEqual(cap, maxCap)) {
          needAdjust = true;
          // clang-format off
          debugPrint(logger_, CTS, "dummy load", 1, "{} => {} "
                     "cap:{:0.2e} idealCap:{:0.2e} delCap:{:0.2e}",
                     driver->getName(), inst->getName(), cap, maxCap,
                     maxCap-cap);
          // clang-format on
        }
      }
    });
  });

  return needAdjust;
}

// Find clock buffers and inverters to use as dummy loads
void TritonCTS::findCandidateDummyCells(
    std::vector<sta::LibertyCell*>& dummyCandidates)
{
  // Add existing buffer list
  for (const std::string& buffer : options_->getBufferList()) {
    odb::dbMaster* master = db_->findMaster(buffer.c_str());

    if (master) {
      sta::Cell* masterCell = network_->dbToSta(master);
      if (masterCell) {
        sta::LibertyCell* libCell = network_->libertyCell(masterCell);
        if (libCell) {
          dummyCandidates.emplace_back(libCell);
        }
      }
    }
  }

  // Add additional inverter cells
  // first, look for inverters with "is_clock_cell: true" cell attribute
  std::vector<sta::LibertyCell*> inverters;
  sta::LibertyLibraryIterator* lib_iter = network_->libertyLibraryIterator();
  while (lib_iter->hasNext()) {
    sta::LibertyLibrary* lib = lib_iter->next();
    for (sta::LibertyCell* inv : *lib->inverters()) {
      if (inv->isClockCell() && resizer_->isClockCellCandidate(inv)) {
        inverters.emplace_back(inv);
        dummyCandidates.emplace_back(inv);
      }
    }
  }
  delete lib_iter;

  // second, look for all inverters with name CLKINV or clkinv
  if (inverters.empty()) {
    sta::PatternMatch patternClkInv("*CLKINV*",
                                    /* is_regexp */ true,
                                    /* nocase */ true,
                                    /* Tcl_interp* */ nullptr);
    lib_iter = network_->libertyLibraryIterator();
    while (lib_iter->hasNext()) {
      sta::LibertyLibrary* lib = lib_iter->next();
      for (sta::LibertyCell* inv :
           lib->findLibertyCellsMatching(&patternClkInv)) {
        if (inv->isInverter() && resizer_->isClockCellCandidate(inv)) {
          inverters.emplace_back(inv);
          dummyCandidates.emplace_back(inv);
        }
      }
    }
    delete lib_iter;
  }

  // third, look for all inverters with name INV or inv
  if (inverters.empty()) {
    sta::PatternMatch patternInv("*INV*",
                                 /* is_regexp */ true,
                                 /* nocase */ true,
                                 /* Tcl_interp* */ nullptr);
    lib_iter = network_->libertyLibraryIterator();
    while (lib_iter->hasNext()) {
      sta::LibertyLibrary* lib = lib_iter->next();
      for (sta::LibertyCell* inv : lib->findLibertyCellsMatching(&patternInv)) {
        if (inv->isInverter() && resizer_->isClockCellCandidate(inv)) {
          inverters.emplace_back(inv);
          dummyCandidates.emplace_back(inv);
        }
      }
    }
    delete lib_iter;
  }

  // abandon attributes & name patterns, just look for all inverters
  if (inverters.empty()) {
    lib_iter = network_->libertyLibraryIterator();
    while (lib_iter->hasNext()) {
      sta::LibertyLibrary* lib = lib_iter->next();
      for (sta::LibertyCell* inv : *lib->inverters()) {
        if (resizer_->isClockCellCandidate(inv)) {
          inverters.emplace_back(inv);
          dummyCandidates.emplace_back(inv);
        }
      }
    }
    delete lib_iter;
  }

  // Sort cells in ascending order of input cap
  std::ranges::sort(
      dummyCandidates,
      [](const sta::LibertyCell* cell1, const sta::LibertyCell* cell2) {
        return (getInputCap(cell1) < getInputCap(cell2));
      });

  if (logger_->debugCheck(utl::CTS, "dummy load", 1)) {
    for (const sta::LibertyCell* libCell : dummyCandidates) {
      // clang-format off
      logger_->debug(CTS, "dummy load",
                     "  {} is a dummy cell candidate with input cap={:0.3e}",
                     libCell->name(), getInputCap(libCell));
      // clang-format on
    }
  }
}

odb::dbInst* TritonCTS::insertDummyCell(
    Clock& clockNet,
    ClockInst* inst,
    const std::vector<sta::LibertyCell*>& dummyCandidates)
{
  ClockSubNet* subNet = driver2subnet_[inst];
  if (subNet->getNumSinks() == options_->getMaxFanout()) {
    return nullptr;
  }
  float deltaCap = inst->getIdealOutputCap() - inst->getOutputCap();
  sta::LibertyCell* dummyCell = findBestDummyCell(dummyCandidates, deltaCap);
  // clang-format off
  debugPrint(logger_, CTS, "dummy load", 1, "insertDummyCell {} at {}",
             inst->getName(), dummyCell->name());
  // clang-format on
  odb::dbInst* dummyInst = nullptr;
  ClockInst& dummyClock = placeDummyCell(clockNet, inst, dummyCell, dummyInst);
  if (driver2subnet_.find(inst) == driver2subnet_.end()) {
    logger_->error(
        CTS, 120, "Subnet was not found for clock buffer {}.", inst->getName());
    return nullptr;
  }
  connectDummyCell(inst, dummyInst, *subNet, dummyClock);
  return dummyInst;
}

ClockInst& TritonCTS::placeDummyCell(Clock& clockNet,
                                     const ClockInst* inst,
                                     const sta::LibertyCell* dummyCell,
                                     odb::dbInst*& dummyInst)
{
  odb::dbMaster* master = network_->staToDb(dummyCell);
  if (master == nullptr) {
    logger_->error(CTS,
                   118,
                   "No phyiscal master cell found for dummy cell {}.",
                   dummyCell->name());
  }
  std::string cellName
      = options_->getDummyLoadPrefix() + std::to_string(dummyLoadIndex_++);
  dummyInst = odb::dbInst::create(block_, master, cellName.c_str());
  dummyInst->setSourceType(odb::dbSourceType::TIMING);
  // V59: Offset dummy load by 2× driver buffer width to prevent physical overlap.
  // 1× bufWidth places dummy at driver edge → DPL can push it 1 site back into
  // driver body. 2× bufWidth gives sufficient margin for DPL micro-adjustment.
  // Dummy load is cap-balance only → position has no timing impact.
  int dummyX = inst->getX();
  int dummyY = inst->getY();
  odb::dbInst* driverDbInst = inst->getDbInst();
  if (driverDbInst && driverDbInst->getMaster()) {
    dummyX += 2 * driverDbInst->getMaster()->getWidth();
  }
  dummyInst->setLocation(dummyX, dummyY);
  dummyInst->setPlacementStatus(odb::dbPlacementStatus::PLACED);
  ClockInst& dummyClock = clockNet.addClockBuffer(
      cellName, master->getName(), dummyX, dummyY);
  // clang-format off
  debugPrint(logger_, CTS, "dummy load", 1, "  placed dummy instance {} at {}",
             dummyInst->getName(), dummyInst->getLocation());
  return dummyClock;
  // clang-format on
}

void TritonCTS::connectDummyCell(const ClockInst* inst,
                                 odb::dbInst* dummyInst,
                                 ClockSubNet& subNet,
                                 ClockInst& dummyClock)
{
  odb::dbInst* sinkInst = inst->getDbInst();
  if (sinkInst == nullptr) {
    logger_->error(
        CTS, 119, "Phyiscal instance {} is not found.", inst->getName());
  }
  odb::dbITerm* iTerm = sinkInst->getFirstOutput();
  odb::dbNet* sinkNet = iTerm->getNet();
  odb::dbITerm* dummyInputPin = getFirstInput(dummyInst);
  dummyInputPin->connect(sinkNet);
  dummyClock.setInputPinObj(dummyInputPin);
  subNet.addInst(dummyClock);
}

void TritonCTS::printClockNetwork(const Clock& clockNet) const
{
  clockNet.forEachSubNet([&](const ClockSubNet& subNet) {
    ClockInst* driver = subNet.getDriver();
    logger_->report("{} has {} sinks", driver->getName(), subNet.getNumSinks());
    subNet.forEachSink([&](const ClockInst* inst) {
      logger_->report("{} -> {}", driver->getName(), inst->getName());
    });
  });
}

void TritonCTS::setAllClocksPropagated()
{
  sta::Sdc* sdc = openSta_->sdc();
  for (sta::Clock* clk : *sdc->clocks()) {
    openSta_->setPropagatedClock(clk);
  }
  estimate_parasitics_->estimateParasitics(est::ParasiticsSrc::placement);
}

void TritonCTS::repairClockNets()
{
  double max_wire_length
      = resizer_->findMaxWireLength(/* don't issue error */ false);
  if (max_wire_length > 0.0) {
    resizer_->repairClkNets(max_wire_length);
  }
}

// Balance macro cell latencies with register latencies.
// This is needed only if special insertion delay handling
// is invoked.
void TritonCTS::balanceMacroRegisterLatencies()
{
  if (!options_->insertionDelayEnabled()) {
    return;
  }

  // Visit builders from bottom up such that latencies are adjusted near bottom
  // trees first
  int totalDelayBuff = 0;

  sta::Corner* corner = openSta_->cmdCorner();
  // convert from per meter to per dbu
  double capPerDBU = estimate_parasitics_->wireClkCapacitance(corner) * 1e-6
                     / block_->getDbUnitsPerMicron();

  for (auto& builder : std::ranges::reverse_view(builders_)) {
    if (builder->getParent() == nullptr && !builder->getChildren().empty()) {
      est::IncrementalParasiticsGuard parasitics_guard(estimate_parasitics_);
      LatencyBalancer balancer = LatencyBalancer(builder.get(),
                                                 options_,
                                                 logger_,
                                                 db_,
                                                 network_,
                                                 openSta_,
                                                 techChar_->getLengthUnit(),
                                                 capPerDBU);
      totalDelayBuff += balancer.run();
    }
  }
  if (totalDelayBuff) {
    logger_->info(CTS, 37, "Total number of delay buffers: {}", totalDelayBuff);
  }
}

void TritonCTS::extractFFGraphFromVerilog(const std::string& verilog_file,
                                          const std::string& output_file)
{
  logger_->report("==========================================");
  logger_->report("Verilog-based FF-to-FF Extraction");
  logger_->report("==========================================");
  logger_->report("Verilog file: {}", verilog_file);
  logger_->report("Output file:  {}", output_file);
  logger_->report("==========================================");

  VerilogFFExtractor extractor(verilog_file, getBlock(), openSta_, network_, logger_);
  // Inject 3D database for tier queries
  if (cts3dDb_) {
    extractor.setCts3DDatabase(cts3dDb_.get());
  }

  // Parse verilog
  extractor.parseVerilog();

  // Extract FF-to-FF edges
  auto edges = extractor.extractFFEdges();

  // Fill timing info from STA (if available)
  // BATCH_TIMING=1: capture-side batched findPathEnds (2-5x faster, exact)
  const char* batch_env = std::getenv("BATCH_TIMING");
  if (batch_env && std::atoi(batch_env) != 0) {
    extractor.fillTimingInfoBatch(edges);
  } else {
    extractor.fillTimingInfo(edges);
  }

  // Write to CSV
  extractor.writeCSV(edges, output_file);

  logger_->report("==========================================");
  logger_->report("Verilog-based extraction complete!");
  logger_->report("Output written to: {}", output_file);
  logger_->report("==========================================");
}

// Exhaustive IO timing edge extraction via STA.
// Replaces Tcl extract_io_timing_edges_phase1a (worst-N) with per-port C++.
void TritonCTS::extractIOTimingEdges(const std::string& verilog_file,
                                     const std::string& output_file)
{
  logger_->report("==========================================");
  logger_->report("3D-CTS IO: Exhaustive IO Timing Edge Extraction");
  logger_->report("==========================================");
  logger_->report("Verilog file: {}", verilog_file);
  logger_->report("Output file:  {}", output_file);

  VerilogFFExtractor extractor(verilog_file, getBlock(), openSta_, network_, logger_);
  if (cts3dDb_) {
    extractor.setCts3DDatabase(cts3dDb_.get());
  }

  // Parse verilog for register name mapping (FF→FF CSV name consistency)
  extractor.parseVerilog();

  // Extract IO edges via STA per-port findPathEnds
  auto io_edges = extractor.extractIOEdges();

  // Write to CSV (same format as Tcl version)
  extractor.writeIOCSV(io_edges, output_file);

  logger_->report("==========================================");
  logger_->report("3D-CTS IO: IO extraction complete!");
  logger_->report("Output written to: {}", output_file);
  logger_->report("==========================================");
}

// BUF_MACRO: ODB+STA-based FF timing graph extraction.
// No Verilog file needed — registers discovered from ODB (Liberty hasSequentials),
// edges extracted via STA findPathEnds. Sees through macro blackboxes.
void TritonCTS::extractFFGraphFromODB(const std::string& output_file)
{
  logger_->report("==========================================");
  logger_->report("ODB+STA-based FF-to-FF Extraction (MACRO)");
  logger_->report("==========================================");
  logger_->report("Output file:  {}", output_file);
  logger_->report("==========================================");

  // Use empty verilog_path — ODB mode doesn't need it
  VerilogFFExtractor extractor("", getBlock(), openSta_, network_, logger_);
  if (cts3dDb_) {
    extractor.setCts3DDatabase(cts3dDb_.get());
  }

  // Step 1: Collect registers from ODB (Liberty hasSequentials)
  extractor.collectRegistersFromODB();

  // Step 2: Extract FF→FF edges via STA findPathEnds
  auto edges = extractor.extractFFEdgesViaSTA();

  // Step 3: Write to CSV (same format as Verilog-based)
  extractor.writeCSV(edges, output_file);

  logger_->report("==========================================");
  logger_->report("ODB+STA-based extraction complete!");
  logger_->report("Output written to: {}", output_file);
  logger_->report("==========================================");
}

// BUF_MACRO: ODB-based IO timing edge extraction.
// No Verilog file needed — uses collectRegistersFromODB() for FF name mapping,
// then calls extractIOEdges() which now falls back to odb_registers_.
void TritonCTS::extractIOTimingEdgesFromODB(const std::string& output_file)
{
  logger_->report("==========================================");
  logger_->report("BUF_MACRO: ODB-based IO Timing Edge Extraction");
  logger_->report("==========================================");
  logger_->report("Output file:  {}", output_file);

  // Use empty verilog_path — ODB mode doesn't need it
  VerilogFFExtractor extractor("", getBlock(), openSta_, network_, logger_);
  if (cts3dDb_) {
    extractor.setCts3DDatabase(cts3dDb_.get());
  }

  // Collect registers from ODB (populates odb_registers_)
  extractor.collectRegistersFromODB();

  // extractIOEdges() now uses odb_registers_ when registers_ is empty
  auto io_edges = extractor.extractIOEdges();

  // Write to CSV (same format as Tcl/Verilog version)
  extractor.writeIOCSV(io_edges, output_file);

  logger_->report("==========================================");
  logger_->report("BUF_MACRO: IO extraction complete!");
  logger_->report("Output written to: {}", output_file);
  logger_->report("==========================================");
}

// ─────────────────────────────────────────────────────────────────────────────
// Pre-CTS skew targets for SG-CTS
// ─────────────────────────────────────────────────────────────────────────────
void TritonCTS::loadSkewTargets(const std::string& csv_path)
{
  // Ensure cts3dDb_ exists (may be called before runTritonCts)
  if (!cts3dDb_) {
    cts3dDb_ = std::make_unique<Cts3DDatabase>(db_, logger_);
    cts3dDb_->populate();
  }
  cts3dDb_->loadSkewTargets(csv_path);
}

// C++ LP solver for skew targeting.
// Replaces Python cts_skew_lp.py — extracts timing graph from STA in-memory,
// solves LP-TNS with OR-Tools GLOP, stores results in skewTargetMap_.
void TritonCTS::solveSkewLp(const std::string& verilog_file,
                            float sigma_local, float sigma_pi,
                            float lambda_reg, float hold_margin,
                            float gamma_wns, float weight_io,
                            bool hard_pi_hold, float max_skew)
{
  // Ensure cts3dDb_ exists
  if (!cts3dDb_) {
    cts3dDb_ = std::make_unique<Cts3DDatabase>(db_, logger_);
    cts3dDb_->populate();
  }

  auto* block = db_->getChip()->getBlock();

  // Create LP solver
  CtsSkewLpSolver lp_solver(block, openSta_, network_, logger_, cts3dDb_.get());

  // Set parameters
  LpParams params;
  params.sigma_local = sigma_local;
  params.sigma_pi = sigma_pi;
  params.lambda_reg = lambda_reg;
  params.hold_margin = hold_margin;
  params.gamma_wns = gamma_wns;
  params.weight_io = weight_io;
  params.hard_pi_hold = hard_pi_hold;
  params.max_skew = max_skew;
  // 2-phase LP params read from env
  // PRE_CTS_ENABLE_2PHASE=1: Phase A minimize W, Phase B minimize Σ V_j
  // PRE_CTS_WNS_ALPHA=0.05: Phase B bound W ≤ W_star*(1+alpha)
  const char* enable_2phase_env = std::getenv("PRE_CTS_ENABLE_2PHASE");
  params.enable_2phase = enable_2phase_env && std::atoi(enable_2phase_env) != 0;
  const char* wns_alpha_env = std::getenv("PRE_CTS_WNS_ALPHA");
  params.wns_alpha = wns_alpha_env ? std::stof(wns_alpha_env) : 0.05f;

  // CTS_DELIVERY_MAX_NS removed — no LP delivery clamp.
  // params.delivery_max stays 0.0 (disabled).

  // V58: Endpoint-TopK pruning for setup edges — controls LP model size.
  // Per capture endpoint j, keep only worst K incoming setup edges.
  // ariane133: 7.9M setup edges → 25K endpoints × 16 = 400K (20× reduction).
  // 0 = disabled (keep all).  Default 16.
  const char* topk_env = std::getenv("CTS_LP_ENDPOINT_TOPK");
  if (topk_env) {
    params.endpoint_topk = std::atoi(topk_env);
  }

  lp_solver.setParams(params);

  // Read Phase 1 CSVs (propagated-clock) instead of live STA (ideal-clock).
  // Phase 1 sub-OpenROAD generates these CSVs with correct propagated-clock slacks.
  // The main session is still at ideal-clock state → live STA gives wrong slacks.
  const char* results_dir = std::getenv("CTS_RESULTS_DIR");
  std::string res_base = results_dir ? std::string(results_dir)
                                     : "./results";
  // Look for the CSV path pattern used by cts_3d.tcl
  const char* graph_csv_env = std::getenv("CTS_FF_TIMING_GRAPH_CSV");
  const char* io_csv_env = std::getenv("CTS_IO_TIMING_EDGES_CSV");

  std::string graph_csv = graph_csv_env ? std::string(graph_csv_env)
                                        : res_base + "/pre_cts_ff_timing_graph.csv";
  // Fallback IO CSV path must include platform/design.
  // Without it, C++ looks at ./results/pre_cts_io_timing_edges.csv (doesn't exist).
  std::string io_csv;
  if (io_csv_env) {
    io_csv = std::string(io_csv_env);
  } else {
    const char* plat = std::getenv("PLATFORM");
    const char* design = std::getenv("DESIGN_NICKNAME");
    if (plat && design) {
      io_csv = res_base + "/" + std::string(plat) + "/" + std::string(design)
               + "/openroad/pre_cts_io_timing_edges.csv";
    } else {
      io_csv = res_base + "/pre_cts_io_timing_edges.csv";
    }
  }

  // Per-clock-net LP removed (V57_CLK, 2026-04-03).
  // Reason: TritonCTS forkRegisterClockNetwork splits a single SDC clock into
  // multiple internal clock nets (e.g., clk, clk_i, clk_i_regs). Per-clock LP
  // treated these as independent clock domains, losing inter-net timing edges.
  // IBEX: single SDC clock → 3 CTS nets → per-clock LP lost 1M+ cross-net edges
  // → each sub-LP still INFEASIBLE (intra-net cycles) + worse timing (info loss).
  // Unified LP is correct: all FFs on the same SDC clock are in one LP.

  // Unified LP path
  lp_solver.readTimingGraphCSV(graph_csv);
  lp_solver.readIOEdgesCSV(io_csv);
  lp_solver.computeBounds();
  lp_solver.applyPiHoldClip();

  LpResult result = lp_solver.solveTns();

  cts3dDb_->clearSkewTargets();
  for (const auto& [ff_name, arrival_ns] : result.arrivals) {
    cts3dDb_->setSkewTarget(ff_name, arrival_ns);
  }

  logger_->info(utl::CTS, 609,
      "3D-CTS LP: Stored {} skew targets ({} non-zero) in C++ memory",
      result.n_ffs, result.n_nonzero_targets);

  // Write targets CSV for downstream consumers.
  // buffer_sizing_lp.py and clock_layer_assignment_v43.tcl read this CSV.
  // Without it, buffer sizing runs with targets=0 → no skew-aware optimization.
  const char* targets_csv_env = std::getenv("CTS_SKEW_TARGETS_CSV");
  std::string targets_csv;
  if (targets_csv_env) {
    targets_csv = std::string(targets_csv_env);
  } else {
    std::string platform = std::getenv("PLATFORM") ? std::getenv("PLATFORM") : "";
    std::string design = std::getenv("DESIGN_NICKNAME") ? std::getenv("DESIGN_NICKNAME") : "";
    if (!platform.empty() && !design.empty()) {
      targets_csv = "./results/" + platform + "/" + design
                    + "/openroad/pre_cts_skew_targets.csv";
    }
  }

  if (!targets_csv.empty()) {
    lp_solver.writeTargetsCSV(targets_csv);
  }

  // BUG #4 FIX: Auto-update TAP depths after post-CTS LP.
  // When CTS_AUTO_UPDATE_TAP=1, the post-CTS LP re-solves with Phase 3
  // propagated-clock timing (correct base latency). We then reconnect FFs
  // to the TAP chain at depths matching the updated LP targets.
  // This fixes the Phase 1↔Phase 3 base latency mismatch:
  //   Phase 1 (balanced): base_latency ~43ps → LP targets based on this
  //   Phase 3 (per-tier):  base_latency ~76ps → LP targets are stale
  //   Post-CTS LP:         base_latency ~76ps → correct targets
  //   → updateTapDepths reconnects FFs to match updated targets.
  const char* auto_tap = std::getenv("CTS_AUTO_UPDATE_TAP");
  if (auto_tap && std::atoi(auto_tap) != 0) {
    // Find targets CSV (just written above)
    std::string tap_csv;
    if (targets_csv_env) {
      tap_csv = std::string(targets_csv_env);
    } else {
      std::string platform = std::getenv("PLATFORM") ? std::getenv("PLATFORM") : "";
      std::string design = std::getenv("DESIGN_NICKNAME") ? std::getenv("DESIGN_NICKNAME") : "";
      if (!platform.empty() && !design.empty()) {
        tap_csv = "./results/" + platform + "/" + design
                  + "/openroad/pre_cts_skew_targets.csv";
      }
    }
    if (!tap_csv.empty()) {
      logger_->info(utl::CTS, 870,
          "Auto-updating TAP depths from post-CTS LP targets: {}", tap_csv);
      updateTapDepths(tap_csv.c_str());
    }
  }
}

// Old C++ buffer sizing LP removed.
// Replaced by BufSizingLpSolver (V53_FM_BUF) with Liberty-based delay estimation.
// See solveBufferSizingLp() below (new signature with CSV paths).

// Estimate per-FF physical achievability bounds for LP-SAFETY.
// Uses ClockLatencyEstimator to compute t_via from HBT parasitic and writes
// bounds CSV consumed by pre_cts_skew_lp.py via --bounds-csv argument.
void TritonCTS::estimateLeafLatencies(const std::string& output_csv,
                                       double max_skew_ns)
{
  // Ensure cts3dDb_ exists and is populated with tier info
  if (!cts3dDb_) {
    cts3dDb_ = std::make_unique<Cts3DDatabase>(db_, logger_);
    cts3dDb_->populate();
  }

  // Load HB via parasitic from tech (set_layer_rc -via hb_layer in Tcl setup)
  // Needed to compute t_via = 0.693 * R_HB * C_HB
  cts3dDb_->loadHybridBondParasitics();

  // Compute per-FF bounds and write CSV
  ClockLatencyEstimator estimator(cts3dDb_.get());
  estimator.estimateAndWrite(output_csv, max_skew_ns);

  const double t_via_ps = estimator.getViaDelayNs() * 1.0e3;
  logger_->report(
      "estimate_leaf_latencies: t_via = {:.2f} ps, max_skew = {:.0f} ps, "
      "written to {}",
      t_via_ps, max_skew_ns * 1.0e3, output_csv);
}

// Update cascaded TAP chain connections in ODB.
// After CTS builds cascaded TAP chains (leaf_buf -> tap1 -> tap2 -> ... -> tap_maxD),
// each FF connects to a specific tap level's output net.  This function re-reads
// LP targets and reconnects FFs to the correct tap depth without rebuilding the tree.
void TritonCTS::updateTapDepths(const char* targets_csv)
{
  // --- 1. Read targets CSV (ff_name, target_ns) ---
  std::ifstream csv(targets_csv);
  if (!csv.is_open()) {
    logger_->warn(CTS, 871, "Cannot open targets CSV: {}", targets_csv);
    return;
  }

  std::unordered_map<std::string, double> targets;
  std::string line;
  std::getline(csv, line);  // skip header
  while (std::getline(csv, line)) {
    auto comma = line.find(',');
    if (comma == std::string::npos)
      continue;
    std::string ff_name = line.substr(0, comma);
    double target_ns = 0.0;
    try {
      target_ns = std::stod(line.substr(comma + 1));
    } catch (...) {
      continue;
    }
    targets[ff_name] = target_ns;
  }
  csv.close();

  // --- 2. Get d_buf from Liberty (same source as cascaded TAP construction) ---
  double d_buf = 0.025;  // fallback
  if (techChar_ && techChar_->getCharBufDelay() > 0) {
    d_buf = techChar_->getCharBufDelay();  // Liberty intrinsic delay (ns)
  }
  // Allow env override for debugging only
  if (const char* e = std::getenv("CTS_PER_FF_BUF_DELAY_NS"))
    d_buf = std::stod(e);

  // --- 3. Get MAX_DEPTH from env var (default 8) ---
  int max_depth = 8;
  if (const char* e = std::getenv("CTS_GROUPED_DELAY_MAX_DEPTH"))
    max_depth = std::atoi(e);

  odb::dbBlock* block = db_->getChip()->getBlock();
  int reconnected = 0, skipped = 0, not_found = 0;

  for (auto& [ff_name, target_ns] : targets) {
    // --- 4a. Compute new depth ---
    int new_depth = std::min(max_depth, std::max(0,
        static_cast<int>(std::round(target_ns / d_buf))));

    // --- 4b. Find FF instance ---
    odb::dbInst* ff_inst = block->findInst(ff_name.c_str());
    if (!ff_inst) {
      not_found++;
      continue;
    }

    // --- 4c. Find CLK input ITerm ---
    odb::dbITerm* clk_iterm = nullptr;
    for (odb::dbITerm* iterm : ff_inst->getITerms()) {
      if (iterm->getSigType() == odb::dbSigType::CLOCK) {
        clk_iterm = iterm;
        break;
      }
    }
    if (!clk_iterm) {
      not_found++;
      continue;
    }

    // --- 4d. Get current net and parse current depth ---
    odb::dbNet* cur_net = clk_iterm->getNet();
    if (!cur_net) {
      not_found++;
      continue;
    }
    std::string cur_net_name = cur_net->getName();

    int cur_depth = 0;
    // Pattern: clknet_ghtap_{branch}_c{idx}_t{depth}_{clk_suffix}
    // e.g. t1_clknet_ghtap_35_c25_t6_clk
    // rfind("_t") finds the _t before depth digits.
    // After depth digits, the clock name suffix follows (e.g. "_clk").
    auto d_pos = cur_net_name.rfind("_t");
    std::string clk_suffix;  // clock name suffix after _t{depth}
    if (d_pos != std::string::npos) {
      try {
        size_t depth_start = d_pos + 2;
        size_t depth_end = depth_start;
        while (depth_end < cur_net_name.size()
               && std::isdigit(cur_net_name[depth_end]))
          depth_end++;
        cur_depth = std::stoi(cur_net_name.substr(depth_start,
                                                   depth_end - depth_start));
        // Preserve suffix after depth digits (e.g. "_clk", "_clk_i_regs")
        clk_suffix = cur_net_name.substr(depth_end);
      } catch (...) {
        cur_depth = 0;
      }
    }
    // If on leaf net (clknet_leaf_*, no _t suffix): depth = 0

    if (new_depth == cur_depth) {
      skipped++;
      continue;
    }

    // --- 4e/f. Find target net and reconnect ---
    std::string target_net_name;

    if (new_depth == 0) {
      // Reconnect to leaf net: trace back through TAP chain to find the leaf net.
      // The depth=1 TAP buffer's input net is the leaf net.
      if (d_pos != std::string::npos) {
        std::string branch_prefix = cur_net_name.substr(0, d_pos);
        // Append clock name suffix to match actual ODB net name
        std::string d1_net_name = branch_prefix + "_t1" + clk_suffix;
        odb::dbNet* d1_net = block->findNet(d1_net_name.c_str());
        if (d1_net) {
          // Find the driving buffer of d1 net, then get its input (= leaf net)
          for (odb::dbITerm* it : d1_net->getITerms()) {
            if (it->getIoType() == odb::dbIoType::OUTPUT) {
              odb::dbInst* tap_buf = it->getInst();
              for (odb::dbITerm* in_it : tap_buf->getITerms()) {
                if (in_it->getIoType() == odb::dbIoType::INPUT
                    && in_it->getSigType() == odb::dbSigType::CLOCK) {
                  odb::dbNet* leaf_net = in_it->getNet();
                  if (leaf_net)
                    target_net_name = leaf_net->getName();
                  break;
                }
              }
              break;
            }
          }
        }
      }
      if (target_net_name.empty()) {
        skipped++;
        continue;
      }
    } else {
      // Reconnect to depth=new_depth TAP net
      if (d_pos != std::string::npos) {
        // Currently on a TAP net — replace _t{old} with _t{new},
        // preserving clock name suffix (e.g. "_clk")
        target_net_name = cur_net_name.substr(0, d_pos)
                          + "_t" + std::to_string(new_depth) + clk_suffix;
      } else {
        // Currently on leaf net (depth=0), find TAP chain driven from this leaf.
        // Scan leaf net's sink ITerms for a TAP buffer instance.
        for (odb::dbITerm* it : cur_net->getITerms()) {
          if (it->getIoType() == odb::dbIoType::INPUT) {
            std::string inst_name = it->getInst()->getName();
            if (inst_name.find("ghtap") != std::string::npos
                || inst_name.find("grp") != std::string::npos) {
              // Found TAP chain entry — get its output net (depth=1)
              for (odb::dbITerm* out_it : it->getInst()->getITerms()) {
                if (out_it->getIoType() == odb::dbIoType::OUTPUT) {
                  odb::dbNet* tap_out = out_it->getNet();
                  if (tap_out) {
                    std::string tap_name = tap_out->getName();
                    auto dp = tap_name.rfind("_t");
                    if (dp != std::string::npos) {
                      // Extract suffix from discovered TAP net
                      size_t ds = dp + 2;
                      while (ds < tap_name.size()
                             && std::isdigit(tap_name[ds]))
                        ds++;
                      std::string tap_suffix = tap_name.substr(ds);
                      target_net_name = tap_name.substr(0, dp)
                                        + "_t" + std::to_string(new_depth)
                                        + tap_suffix;
                    }
                  }
                  break;
                }
              }
              break;
            }
          }
        }
        if (target_net_name.empty()) {
          skipped++;
          continue;
        }
      }
    }

    // Find target net in ODB
    odb::dbNet* target_net = block->findNet(target_net_name.c_str());
    if (!target_net) {
      // Target TAP depth net doesn't exist (chain not deep enough)
      logger_->warn(CTS, 872,
          "TAP net {} not found for FF {} (depth {}->{})",
          target_net_name, ff_name, cur_depth, new_depth);
      not_found++;
      continue;
    }

    // Reconnect: disconnect from old net, connect to new net
    clk_iterm->disconnect();
    clk_iterm->connect(target_net);
    reconnected++;
  }

  logger_->info(CTS, 873,
      "TAP depth update: {} FFs reconnected, {} unchanged, {} not found "
      "(d_buf={:.4f}ns, max_depth={})",
      reconnected, skipped, not_found, d_buf, max_depth);
}

// C++ LP-based buffer sizing with Liberty delays
void TritonCTS::solveBufferSizingLp(const char* output_csv,
                                     const char* skew_targets_csv,
                                     double hold_weight, double reg_weight,
                                     double skew_weight,
                                     double setup_margin_ps, double hold_margin_ps)
{
  BufSizingLpSolver solver(db_->getChip()->getBlock(),
                           openSta_, logger_);

  BufSizingParams params;
  params.hold_weight = hold_weight;
  params.reg_weight = reg_weight;
  params.skew_weight = skew_weight;
  params.setup_margin_ps = setup_margin_ps;
  params.hold_margin_ps = hold_margin_ps;
  solver.setParams(params);

  solver.extractBuffers();
  solver.measureLibraryDelays();
  solver.extractTimingEdges();

  if (skew_targets_csv && std::string(skew_targets_csv) != "") {
    solver.loadSkewTargets(skew_targets_csv);
  }

  solver.solve();
  solver.writeResultsCSV(output_csv);
}

}  // namespace cts
