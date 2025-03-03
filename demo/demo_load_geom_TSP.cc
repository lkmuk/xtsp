#include "xtsp/core/complete_graph.h"

#include <spdlog/spdlog.h>
#include <iostream>

int main(int argc, char const *argv[])
{
  spdlog::set_level(spdlog::level::debug);
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
    SPDLOG_INFO("result: M = {:d}", g.numClusters());
  }
  SPDLOG_INFO("result: dimension of each XY = {:d}", g.nDim());
  SPDLOG_INFO("result: XY data shown below");

  std::cout << g.getXy();
  std::cout << std::endl;

  return 0;
}
