#include "xtsp/core/tour.h"

#include <gtest/gtest.h>
#define SPDLOG_ACTIVE_LEVEL SPDLOG_INFO
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>


TEST(hamiltonianPermTour, successfulInit)
{
  xtsp::PermTour tour({3,2,1,0,4});
  EXPECT_EQ(tour.size(), 5);
}

TEST(hamiltonianTour, directAccessSuccess)
{
  xtsp::PermTour tour({3,2,1,0,4});
  EXPECT_EQ(tour.getVertex_(0), 3);
  EXPECT_EQ(tour.getVertex_(1), 2);
  EXPECT_EQ(tour.getVertex_(2), 1);
  EXPECT_EQ(tour.getVertex_(3), 0);
  EXPECT_EQ(tour.getVertex_(4), 4);
}

TEST(hamiltonianTour, circularBufAccessSuccess)
{
  xtsp::PermTour tour({3,2,1,0,4});
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
  xtsp::PermTour tour({2, 3, 0, 1});

  const int expectedCost = 1000+10+5+99999;
  EXPECT_EQ(evalTour(tour, g), expectedCost);
}

// suffices to test the assertion is really execution.
// We don't test other faulty calls (i.e., mismatchSize/missing/too much, erroneous values)
TEST(hamiltonianTour, catchDuplicate)
{
  const std::vector<size_t> testArr{3, 2, 4, 0, 2};

  // change this
  const std::string expectedMsg =
      "Invalid tour because city 2 appears at least twice";
  try
  {
    // change this
    xtsp::PermTour tour (testArr, 5);
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

  auto gtour = xtsp::GeneralizedTour::fromPermutation({3,2,4}, clustering);

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

  auto gtour = xtsp::GeneralizedTour::fromPermutation({3,5,4}, clustering);
  EXPECT_EQ(gtour.getTour()->getVertex(1+9), 5);
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
    auto gtour = xtsp::GeneralizedTour::fromPermutation (testPerm, clustering);
    FAIL() << fmt::format("Should have thrown an exception: {}", expectedMsg);
  }
  catch (const std::invalid_argument &actualException)
  {
    EXPECT_EQ(actualException.what(), expectedMsg);
  }
}


class SimpleTwoEdgeExchange : public testing::Test
{
protected:
  // only FYI not really critical in this context
  //
  // In fact, this is the simplest uncrossing example
  // with 4 points on a unit circle
  Eigen::Matrix<float, 4, 2> xy {
    {1, 0},
    {0, 1},
    {-1, 0},
    {0, -1}
  };

  xtsp::ImplicitCompleteGraph<float> g = xtsp::ImplicitCompleteGraph<float>(xy);

  const float expectedOldTourCost = std::sqrt(2)*2 + 2*2;
  const float expectedNewTourCost = std::sqrt(2)*4;
  std::vector<size_t> originalTourDat = {2, 0, 3, 1};
  xtsp::PermTour tour =  xtsp::PermTour(originalTourDat);
  // by eyeballing, we identify a Two-opt
  const size_t rankA = 1; // vertex 0
  // const size_t rankB = 2; // vertex 3 (only FYI)
  // ensure rankC >= rank A + 1 (otherwise += numVertices)
  const size_t rankC = 3; // vertex 1 (in general the rankC >= rankB)
  // const size_t rankD = 0; // vertex 2 (only FYI)

  const size_t vA = 0;
  const size_t vC = 1;

  void SetUp() override
  {
    SPDLOG_DEBUG("Tour sequence (before): {:d}-{:d}-{:d}-{:d}-",
      tour.getVertex_(0), 
      tour.getVertex_(1), 
      tour.getVertex_(2), 
      tour.getVertex_(3)
    );
    EXPECT_FLOAT_EQ(xtsp::evalTour(tour, g), expectedOldTourCost) 
      << "bug in the test implementation";
  }

  void checkResultAfterTwoOpt() const
  {
    tour.isHamiltonian();

    EXPECT_FLOAT_EQ(xtsp::evalTour(tour, g), expectedNewTourCost); // check if it's legit
    SPDLOG_DEBUG("Tour sequence (after):  {:d}-{:d}-{:d}-{:d}-",
        tour.getVertex_(0), 
        tour.getVertex_(1), 
        tour.getVertex_(2), 
        tour.getVertex_(3)
    );
  }
};

TEST_F(SimpleTwoEdgeExchange, permTourRankBased)
{
  spdlog::set_level(spdlog::level::debug);
  tour.exchangeTwoEdges_rankBased(rankA, rankC);
  checkResultAfterTwoOpt();
}

TEST_F(SimpleTwoEdgeExchange, permTourNative)
{
  spdlog::set_level(spdlog::level::debug);
  tour.exchangeTwoEdges(vA, vC);
  checkResultAfterTwoOpt();
}
