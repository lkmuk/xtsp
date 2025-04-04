#include "xtsp/core/tour_alternatives.h"
#include "xtsp/initialization/insertion.h"
#include "xtsp/core/utils.h"
#include "xtsp/local_search/kopt.h"

#include <gtest/gtest.h>
#include <filesystem>
#include <map>
#include <chrono>

#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>

static const auto dataDir = std::filesystem::path(
  __FILE__).parent_path().parent_path()/"tests"/"dataset";

int main(int argc, char const *argv[])
{
  spdlog::set_level(spdlog::level::info);
  if (argc != 2)
  {
    SPDLOG_ERROR(
      "Wrong usage, please specify the file of a " 
      "geometry TSP instance in TSPLIB format.");
    return -1;
  }

  auto g = xtsp::ImplicitCompleteGraph<float>::loadFromTsplibFile(argv[1]);
  SPDLOG_INFO("Successfully loaded a Geometric TSP instance from a file.");
  SPDLOG_INFO("result: N = {:d}", g.numVertices());
  /// @todo assert is_symmetric
  if (g.isClustered())
  {
    SPDLOG_WARN("ignoring the clustering");
  }
  float upscale = 1;
  auto gExplicit = g.explicitize(upscale);

  SPDLOG_INFO("Begin tour construction");
  auto t0 = std::chrono::high_resolution_clock::now();

  auto tour = xtsp::algo::farthestInsertion(gExplicit, 0);
  // std::mt19937 rng {324};// (std::random_device{}());
  // std::vector<size_t> initPerm;
  // xtsp::utils::genPermutation(rng, gExplicit.numVertices(), initPerm);
  // auto tour = xtsp::PermTour(initPerm);

  auto t1 = std::chrono::high_resolution_clock::now();
  if (!tour.isHamiltonian())
    SPDLOG_ERROR("Resultant tour is not Hamiltonian!");
  auto costInit = xtsp::evalTour(tour, gExplicit)/upscale;
  SPDLOG_INFO(
    "Finished tour construction: cost = {}, + {:.3f} s", 
    costInit, (t1-t0).count()*1e-9);
  
  auto t2 = std::chrono::high_resolution_clock::now();

  // typically still faster despite the conversion overhead
  // e.g., in pr144
  // auto tourRefined = tour;
  xtsp::AdjTabTour tourRefined (tour.getSeqMutableRef__()); 
  xtsp::algo::PriorityTwoOptFinder solver(tourRefined, gExplicit);
  auto twoOptOutcome = solver.solve(tourRefined, gExplicit, 100, true);

  auto t3 = std::chrono::high_resolution_clock::now();
  auto costNew = xtsp::evalTour(tourRefined, gExplicit)/upscale;
  {
    auto gap = costNew + twoOptOutcome.improvement()/upscale -costInit;
    if (std::abs(gap) > 1/upscale)
      SPDLOG_WARN("Inconsistent cost value, possibly bug in the two-opt algorithm, or maybe overflow in the initial tour (especially for randomly generated tours)");
  }
  SPDLOG_INFO(
    "Finished two-opt           : cost = {}, + {:.3f} s",
    costNew, (t3-t2).count()*1e-9);
  if (twoOptOutcome.confirmedTwoOpt())
    SPDLOG_INFO("Final tour is two-opt.");

  auto probName = std::filesystem::path(argv[1]).stem();
  std::string fname = fmt::format("{}.{}.tour", 
    probName.string(), costNew);
  SPDLOG_INFO("Writing the final tour to {}", fname);
  tourRefined.saveTsplib(fname, probName);
  return 0;
}
