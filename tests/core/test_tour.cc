#include "xtsp/core/tour.h"
#include "xtsp/core/tour_alternatives.h"

#include <gtest/gtest.h>
#define SPDLOG_ACTIVE_LEVEL SPDLOG_INFO
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>


TEST(hamiltonianPermTour, successfulInit)
{
  xtsp::PermTour tour({3,2,1,0,4});
  EXPECT_EQ(tour.size(), 5);
}

TEST(hamiltonianAdjTabTour, successfulInit)
{
  xtsp::AdjTabTour tour({3,2,1,0,4});
  EXPECT_EQ(tour.size(), 5);
  EXPECT_EQ(tour.getDepotId(), 3);

  EXPECT_EQ(tour.next(3), 2);
  EXPECT_EQ(tour.next(2), 1);
  EXPECT_EQ(tour.next(1), 0);
  EXPECT_EQ(tour.next(0), 4);
  EXPECT_EQ(tour.next(4), 3); // across the cut/depot
  EXPECT_TRUE(tour.isOneStepAhead(2, 1));
  EXPECT_TRUE(tour.isOneStepAhead(4, 3)); // across the cut/depot

  EXPECT_TRUE(tour.allElementsAreValid());
  EXPECT_TRUE(tour.isHamiltonian());
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
  xtsp::PermTour tourPermType =  xtsp::PermTour(originalTourDat);
  xtsp::AdjTabTour tourAdjType = xtsp::AdjTabTour(originalTourDat);
  // by eyeballing, we identify a Two-opt
  const size_t rankA = 2;
  const size_t rankC = 0+g.numVertices();
  const size_t vA = 3;
  const size_t vC = 2;

  void SetUp() override
  {
    spdlog::set_level(spdlog::level::debug);
    SPDLOG_DEBUG("Tour sequence (before): {:d}-{:d}-{:d}-{:d}-",
      tourPermType.getVertex_(0), 
      tourPermType.getVertex_(1), 
      tourPermType.getVertex_(2), 
      tourPermType.getVertex_(3)
    );
    EXPECT_FLOAT_EQ(xtsp::evalTour(tourPermType, g), expectedOldTourCost) 
      << "bug in the test implementation";
  }

};

TEST_F(SimpleTwoEdgeExchange, permTourRankBased)
{
  tourPermType.exchangeTwoEdges_rankBased(rankA, rankC, true);
  EXPECT_TRUE(tourPermType.isHamiltonian());
  EXPECT_FLOAT_EQ(xtsp::evalTour(tourPermType, g), expectedNewTourCost);
  EXPECT_EQ(tourPermType.getVertex_(0), 1);
  EXPECT_EQ(tourPermType.getVertex_(1), 0);
  EXPECT_EQ(tourPermType.getVertex_(2), 3);
  EXPECT_EQ(tourPermType.getVertex_(3), 2);
  
}

TEST_F(SimpleTwoEdgeExchange, permTour)
{
  tourPermType.exchangeTwoEdges(vA, vC, true);
  EXPECT_TRUE(tourPermType.isHamiltonian());
  EXPECT_FLOAT_EQ(xtsp::evalTour(tourPermType, g), expectedNewTourCost);
  EXPECT_EQ(tourPermType.getVertex_(0), 1);
  EXPECT_EQ(tourPermType.getVertex_(1), 0);
  EXPECT_EQ(tourPermType.getVertex_(2), 3);
  EXPECT_EQ(tourPermType.getVertex_(3), 2);
}

TEST_F(SimpleTwoEdgeExchange, adjTabTour)
{

  EXPECT_TRUE(tourAdjType.isHamiltonian());
  EXPECT_FLOAT_EQ(xtsp::evalTour(tourAdjType, g), expectedOldTourCost);
  
  EXPECT_TRUE(tourAdjType.isTwoPlusStepsAhead(vA,vC));
  EXPECT_TRUE(tourAdjType.isTwoPlusStepsAhead(vC,vA));

  // test subject
  tourAdjType.exchangeTwoEdges(vA, vC, true);
  
  EXPECT_TRUE(tourAdjType.isHamiltonian());
  std::cout << "After 2-opt: ";
  tourAdjType.print(std::cout);
  std::cout << std::endl;
  // possible perm:  2 - 1 - 0 - 3 
  //      1 - 0 - 3 - 2
  EXPECT_EQ(tourAdjType.next(0), 3);
  EXPECT_EQ(tourAdjType.next(1), 0);
  EXPECT_EQ(tourAdjType.next(2), 1);
  EXPECT_EQ(tourAdjType.next(3), 2);
  EXPECT_FLOAT_EQ(xtsp::evalTour(tourAdjType, g), expectedNewTourCost);
}
