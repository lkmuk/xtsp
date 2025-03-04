#include "xtsp/initialization/insertion.h"

#include <gtest/gtest.h>
#include <filesystem>
#include <map>

#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>

static const auto dataDir = std::filesystem::path(
  __FILE__).parent_path().parent_path()/"dataset";

TEST(farthestInsertion, pr144)
{
  // http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp95.pdf
  // instance name and true min (using integer arithmetic!!!)
  std::map<std::string, float> oracle; 
  oracle["pr144"] = 58537;
  oracle["u1817"] = 57201;

  for (const auto&[name, trueMin] : oracle)
  {
    // setup
    spdlog::set_level(spdlog::level::warn);
    auto fpath = dataDir/fmt::format("{}.tsp", name);
    auto g = xtsp::ImplicitCompleteGraph<float>::loadFromTsplibFile(fpath.c_str());

    spdlog::set_level(spdlog::level::debug);
    SPDLOG_INFO("begin construction for {}", name); // see the timestamp
    auto tour = xtsp::algo::farthestInsertion(g);
    SPDLOG_INFO("end construction with tour cost: {:.3f} (+{:.3}% true min)", 
      tour.getCost(), (tour.getCost()-trueMin)/trueMin*100); 
  }
}