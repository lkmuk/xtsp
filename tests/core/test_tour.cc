#include "xtsp/core/tour.h"

#include <gtest/gtest.h>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>


TEST(hamiltonianTour, successfulInit)
{
  xtsp::Tour<int> tour({3,2,1,0,4}, 5);
  EXPECT_EQ(tour.size(), 5);
}

TEST(hamiltonianTour, directAccessSuccess)
{
  xtsp::Tour<int> tour({3,2,1,0,4}, 5);
  EXPECT_EQ(tour.getVertex_(0), 3);
  EXPECT_EQ(tour.getVertex_(1), 2);
  EXPECT_EQ(tour.getVertex_(2), 1);
  EXPECT_EQ(tour.getVertex_(3), 0);
  EXPECT_EQ(tour.getVertex_(4), 4);
}

TEST(hamiltonianTour, circularBufAccessSuccess)
{
  xtsp::Tour<int> tour({3,2,1,0,4}, 5);
  EXPECT_EQ(tour.getVertex(0), 3);
  EXPECT_EQ(tour.getVertex(1+5), 2);
  EXPECT_EQ(tour.getVertex(1+5+5), 2);
  EXPECT_EQ(tour.getVertex(4+5*10), 4);
}

TEST(hamiltonianTour, costEval)
{
  Eigen::Matrix4i costs
  {
    {0, 111, 111, 111},
    {5, 0, 111, 111},
    {1, 99999, 0, 111}, 
    {10, 100, 1000, 0}
  };
  xtsp::CompleteGraph<int> g(true, costs);
  // Note that the tour cost won't be as high as 99999 as long as the edge (1-2) is not in the tour
  xtsp::Tour<int> tour({2, 3, 0, 1}, 4);

  const int expectedCost = 1000+10+5+99999;
  EXPECT_EQ(tour.evalTour(g), expectedCost);

  EXPECT_EQ(tour.getCost(), expectedCost);
}

// suffices to test the assertion is really execution.
// We don't test other faulty calls (i.e., mismatchSize/missing/too much, erroneous values)
TEST(hamiltonianTour, catchDuplicate)
{
  const std::vector<size_t> testArr{3, 2, 4, 0, 2};

  // change this
  const std::string expectedMsg =
      "Invalid tour because city 2 appears at least twice (at tour positions 1 and 4.)";
  try
  {
    // change this
    xtsp::Tour<int> tour (testArr, 5);
    FAIL() << fmt::format("Should have thrown an exception: {}", expectedMsg);
  }
  catch (const std::invalid_argument &actualException)
  {
    EXPECT_EQ(actualException.what(), expectedMsg);
  }
}


TEST(generalizedTour, successfulInit)
{
  xtsp::Clustering clustering_(5, {{2},{0,4},{1,3}});
  auto clustering = std::make_shared<xtsp::Clustering>(clustering_);
  xtsp::GeneralizedTour<int> gtour({3,2,4}, clustering);

  EXPECT_EQ(gtour.numClusters(), clustering->numClusters());
  EXPECT_EQ(gtour.numVertices(), clustering->numVertices());

  EXPECT_EQ(gtour.getClusterIdByRank(0), 2);
  EXPECT_EQ(gtour.getClusterIdByRank(1), 0);
  EXPECT_EQ(gtour.getClusterIdByRank(2), 1);

  // reverse lookup from the internal data
  EXPECT_EQ(gtour.findClusterRankById(0), 1);
  EXPECT_EQ(gtour.findClusterRankById(1), 2);
  EXPECT_EQ(gtour.findClusterRankById(2), 0);


  EXPECT_EQ(gtour.getVertexByClusterId(0), 2);
  EXPECT_EQ(gtour.getVertexByClusterId(1), 4);
  EXPECT_EQ(gtour.getVertexByClusterId(2), 3);
}
TEST(generalizedTour, circularBufAccessSuccess)
{
  xtsp::Clustering clustering_(7, {{0,1,2,4}, {5,6}, {3}});
  auto clustering = std::make_shared<xtsp::Clustering>(clustering_);


  xtsp::GeneralizedTour<int> gtour({3,5,4}, clustering);
  EXPECT_EQ(gtour.getVertex(1+9), 5);
}
TEST(generalizedTour, costValueCanBeUpdated)
{
  xtsp::Clustering clustering_(7, {{0,1,2,4}, {5,6}, {3}});
  auto clustering = std::make_shared<xtsp::Clustering>(clustering_);

  xtsp::GeneralizedTour<int> gtour({3,5,4}, clustering);
  gtour.getCostMutableRef() = 111;
  gtour.getCostMutableRef() += 32;
  EXPECT_EQ(gtour.getCost(), 143);
}

TEST(generalizedTour, catchRevisitClusterWithDifferentVertices)
{
  xtsp::Clustering clustering_(6, {{2, 3, 0}, {4}, {1, 5}});
  auto clustering = std::make_shared<xtsp::Clustering>(clustering_);
  const std::vector<size_t> testPerm {3,5,0};


  // change this
  const std::string expectedMsg =
      "Invalid generalized tour because cluster 0 appears at least twice";
  try
  {
    // change this
    xtsp::GeneralizedTour<int> gtour (testPerm, clustering);
    FAIL() << fmt::format("Should have thrown an exception: {}", expectedMsg);
  }
  catch (const std::invalid_argument &actualException)
  {
    EXPECT_EQ(actualException.what(), expectedMsg);
  }
}
