#include "xtsp/core/clustering.h"

#include <gtest/gtest.h>
// #include "../expect_throw_and_msg.h"
#include <spdlog/spdlog.h>

namespace fmt = spdlog::fmt_lib;


TEST(clustering, successfulInit)
{
  spdlog::set_level(spdlog::level::debug);

  xtsp::Clustering clustering(
      8,
      {
          {6, 4},    // cluster 0
          {2, 1, 5}, // cluster 1
          {0, 3},    // cluster 2
          {7}        // cluster 3
      });

  EXPECT_EQ(clustering.numClusters(), 4);
  EXPECT_EQ(clustering.numVertices(), 8);

  const std::vector<size_t> expectedRes = {
      2, 1, 1, 2, 0, 1, 0, 3};
  ASSERT_EQ(expectedRes.size(), clustering.numVertices())
      << "Bug in the test implementation";
  for (size_t i = 0; i < clustering.numVertices(); ++i)
  {
    EXPECT_EQ(clustering.getClusterId(i), expectedRes[i]) << "at i = " << i;
  }
}

TEST(clustering, queryCluster)
{
  const std::vector<size_t> lastCluster = {2, 3, 1};
  xtsp::Clustering cluster(4, {{0}, lastCluster});

  for (size_t i = 0; i < lastCluster.size(); ++i)
  {
    EXPECT_EQ(cluster.getMembers(1)[i], lastCluster[i]) << "at i = " << i;
  }
}

TEST(clustering, catchInvalidVertex)
{
  // change this
  const std::string expectedMsg = 
    "Invalid clustering because the maximum vertex index of N-1, i.e., 4 is exceeded";
  try
  {
      // change this
      xtsp::Clustering cluster(4, {{0}, {1, 2, 3, 4}});
      FAIL() << fmt::format("Should have thrown an exception: {}", expectedMsg);
  } 
  catch (const std::invalid_argument& actualException)
  {
      EXPECT_EQ(actualException.what(), expectedMsg);
  }

  // expect_throw_thisMsg<>(
  //   [](){;},
  //   "");
}

TEST(clustering, catchMissingVertex)
{
  // change this
  const std::string expectedMsg = 
    "Invalid clustering because vertex 2 is not assigned to any cluster.";
  try
  {
      // change this
      xtsp::Clustering cluster(4, {{0}, {1, 3}});
      FAIL() << fmt::format("Should have thrown an exception: {}", expectedMsg);
  } 
  catch (const std::invalid_argument& actualException)
  {
      EXPECT_EQ(actualException.what(), expectedMsg);
  }
}

TEST(clustering, catchTwoPartitionsOverlap)
{
  /// vertex 2 assigned multiple times

  // change this
  const std::string expectedMsg = 
    "Invalid clustering because each vertex shall be assigned to exactly one cluster "
    "but you assign vertex 2 to both clusters 0 and 2.";
  try
  {
      // change this
      xtsp::Clustering cluster(4, {{0, 2}, {3}, {1, 2}});
      FAIL() << fmt::format("Should have thrown an exception: {}", expectedMsg);
  } 
  catch (const std::invalid_argument& actualException)
  {
      EXPECT_EQ(actualException.what(), expectedMsg);
  }
}