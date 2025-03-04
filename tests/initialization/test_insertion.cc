#include "xtsp/initialization/insertion.h"

#include <gtest/gtest.h>
#include <filesystem>
#include <map>
#include <chrono>

#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>

static const auto dataDir = std::filesystem::path(
  __FILE__).parent_path().parent_path()/"dataset";

class FarthestInsertionTest : public testing::Test
{
protected:
  // http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp95.pdf
  // instance name and true min (using integer arithmetic!!!)
  std::map<std::string, float> oracle; 

  void SetUp () override
  {
    oracle["pr144"] = 58537;
    oracle["u1817"] = 57201;
  }
  std::string filePath (const std::string& probName) const
  {
    return dataDir/fmt::format("{}.tsp", probName);
  }
  float percentGapTrueMin(const std::string& probName, float computedCost)
  {
    return ((computedCost - oracle[probName])/oracle[probName])*100;
  }
};

TEST_F(FarthestInsertionTest, pr114)
{
  std::string name = "pr144";
  // setup
  spdlog::set_level(spdlog::level::warn);
  auto fpath = filePath(name);
  auto g = xtsp::ImplicitCompleteGraph<float>::loadFromTsplibFile(fpath.c_str());

  spdlog::set_level(spdlog::level::debug);
  SPDLOG_INFO("begin construction for {}", name); // see the timestamp
  auto tour = xtsp::algo::farthestInsertion(g);
  SPDLOG_INFO("end construction with tour cost: {:.3f} (+{:.3}% true min)", 
    tour.getCost(), percentGapTrueMin(name, tour.getCost()));

  SPDLOG_INFO("begin explicitization for {} (scale = 1)", name);
  auto gExpanded = g.explicitize(1);
  SPDLOG_INFO("begin construction for {}", name); // see the timestamp
  auto tourInt = xtsp::algo::farthestInsertion(gExpanded);
  SPDLOG_INFO("end construction with tour cost: {:d} (+{:.3f}% true min)", 
    tourInt.getCost(), percentGapTrueMin(name, tourInt.getCost()));

  EXPECT_LT(percentGapTrueMin(name, tourInt.getCost()), 20); // typically won't be worse than 20 for this kind of toy problem

}

TEST_F(FarthestInsertionTest, pr114ImpactFirstPick)
{
  std::string name = "pr144";
  // setup
  spdlog::set_level(spdlog::level::warn);
  auto fpath = filePath(name);
  auto g = xtsp::ImplicitCompleteGraph<float>::loadFromTsplibFile(fpath.c_str());

  spdlog::set_level(spdlog::level::debug);
  ASSERT_EQ(g.numVertices(), 144); // required assumption for what follows
  std::vector<size_t> listPick1 = {1, 5, 10, 15, 20, 67, 80, 100, 100};
  Eigen::Vector<float, -1> tourCost (listPick1.size());
  int cnt = -1;
  for (const size_t vPick1 : listPick1)
  {
    auto tour = xtsp::algo::farthestInsertion(g, vPick1);
    tourCost(++cnt) = tour.getCost(); 
  }
  SPDLOG_INFO("Finished {:d} calls to farthestInsertion for {}", listPick1.size(), name);
  SPDLOG_INFO("Gap of tour cost to true min: {:.3f} (best), {:.3f} (mean), {:.3f} (worst)",
    percentGapTrueMin(name, tourCost.minCoeff()), 
    percentGapTrueMin(name, tourCost.mean()), 
    percentGapTrueMin(name, tourCost.maxCoeff()));  
}

TEST_F(FarthestInsertionTest, u1817RuntimeTest)
{
  std::string name = "u1817";
  // setup
  spdlog::set_level(spdlog::level::warn);
  auto fpath = filePath(name);
  auto g = xtsp::ImplicitCompleteGraph<float>::loadFromTsplibFile(fpath.c_str());
  SPDLOG_INFO("begin explicitization for {} (upscaling: 1)",  name);
  auto gExpanded = g.explicitize(1);

  spdlog::set_level(spdlog::level::debug);
  auto t0 = std::chrono::high_resolution_clock::now();
  SPDLOG_INFO("begin construction for {}", name); // see the timestamp
  auto tourInt = xtsp::algo::farthestInsertion(gExpanded, 0);
  auto tf = std::chrono::high_resolution_clock::now();
  SPDLOG_INFO("end construction with tour cost: {:d} (+{:.3f}% true min) in {:.3f} s", 
    tourInt.getCost(), 
    percentGapTrueMin(name, tourInt.getCost()),
    (tf-t0).count()*1e-9); 
}

