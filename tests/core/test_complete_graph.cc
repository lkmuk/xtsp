#include "xtsp/core/complete_graph.h"

#include <gtest/gtest.h>
#include <spdlog/spdlog.h>
#include <sstream>

class SimpleCompleteGraphInt : public testing::Test
{
protected:
  // which is deliberately asymetric
  Eigen::Matrix<int, 4,4> costMatrix
  {
    {00, 3, 4, 15},
    {15, 0, 2, 1},
    {12, 3, 0, 8},
    {23, 7, 6, 0},
  };
};

TEST_F(SimpleCompleteGraphInt, asymmetricCostQuery)
{
  xtsp::CompleteGraph<int> g (false, this->costMatrix);

  EXPECT_FALSE(g.isSymmetric());
  EXPECT_FALSE(g.isClustered());

  for (size_t i = 0; i < (size_t)costMatrix.rows(); ++i)
    for (size_t j = 0; j < (size_t)costMatrix.cols(); ++j)
      EXPECT_EQ(g.getEdgeCost(i,j), this->costMatrix(i,j)) 
        << "wrong edge cost value at (i,j) = (" << i << ", " << j << ")";
}

TEST_F(SimpleCompleteGraphInt, symmetricCostQuery)
{
  xtsp::CompleteGraph<int> g (true, this->costMatrix);

  EXPECT_TRUE(g.isSymmetric());
  EXPECT_FALSE(g.isClustered());

  for (size_t i = 0; i < (size_t)costMatrix.rows(); ++i)
    for (size_t j = 0; j < (size_t)costMatrix.cols(); ++j)
      EXPECT_EQ(g.getEdgeCost(i,j), g.getEdgeCost(j,i)) 
        << "not symmetric at (i,j) = (" << i << ", " << j << ")";

  // copied the lower-triangular?
  for (size_t i = 0; i < (size_t)costMatrix.rows(); ++i)
    for (size_t j = 0; j < i; ++j)
      EXPECT_EQ(g.getEdgeCost(i,j), this->costMatrix(i,j)) 
        << "wrong edge cost value at (i,j) = (" << i << ", " << j << ")";
}

class VerySmallImplicitCompleteGraph2D : public testing::Test
{
protected:
  const Eigen::Matrix<float,-1, 2> xy
  {
    {-1, 1},
    {2, -3},
    {0, 0},
  };
};

TEST_F(VerySmallImplicitCompleteGraph2D, normTy2)
{
  xtsp::ImplicitCompleteGraph<float> g(xy, nullptr);
  EXPECT_EQ(g.nDim(), 2);
  EXPECT_FLOAT_EQ(g.getEdgeCost(0, 1), 5.);
  EXPECT_FLOAT_EQ(g.getEdgeCost(1, 2), std::sqrt(4+9));
  EXPECT_FLOAT_EQ(g.getEdgeCost(2, 0), std::sqrt(2));
}
TEST_F(VerySmallImplicitCompleteGraph2D, normTy1)
{
  xtsp::ImplicitCompleteGraph<float> g(xy, nullptr, 1);
  EXPECT_FLOAT_EQ(g.getEdgeCost(0, 1), 7.);
  EXPECT_FLOAT_EQ(g.getEdgeCost(1, 2), 5.);
  EXPECT_FLOAT_EQ(g.getEdgeCost(2, 0), 2.);
}
TEST_F(VerySmallImplicitCompleteGraph2D, normTyInf)
{
  xtsp::ImplicitCompleteGraph<float> g(xy, nullptr, 0);
  EXPECT_FLOAT_EQ(g.getEdgeCost(0, 1), 4.);
  EXPECT_FLOAT_EQ(g.getEdgeCost(1, 2), 3.);
  EXPECT_FLOAT_EQ(g.getEdgeCost(2, 0), 1.);
}


class SimpleImplicitCompleteGraph3D : public testing::Test
{
protected:
  const Eigen::Matrix<float,-1, 3> xyz
  {
    {0.0, 0.2, 0.4},
    {-0.1, 0.3, 0.2},
    {-30.1, 34.1, 54.3},
    {4356.3, 0.2, -20.0},
  };
  const size_t numPoints = xyz.rows();
};
TEST_F(SimpleImplicitCompleteGraph3D, costQuery)
{
  xtsp::ImplicitCompleteGraph<float> g (this->xyz);
  EXPECT_EQ(g.nDim(), 3);
  EXPECT_TRUE(g.isSymmetric());
  EXPECT_FALSE(g.isClustered());

  // To be sure the expectedVal is correct 
  // (1) switch to debug level (2) rebuild the test 
  //
  spdlog::set_level(spdlog::level::info);
  
  std::stringstream ss;
  ss << g.getXy();
  spdlog::debug("copied XYZ value \n{}", ss.str());


  for (size_t i = 0; i < numPoints; ++i)
  {
    for (size_t j = 0; j < numPoints; ++j)
    {
      float expectedVal = (xyz.row(i)-xyz.row(j)).norm();
      spdlog::debug("dist(i={:d}, j={:d}) = {:f}", i,j, expectedVal);
      EXPECT_EQ(g.getEdgeCost(i,j), expectedVal)
        << "wrong edge cost at (i,j) = (" << i << ", " << j << ")";
    }
  }
}

TEST(clusterSuperGraph, usingAveraging)
{
  // setup the test
  Eigen::Matrix<float, 9, 2> xyFull;
  xyFull <<
    1,3,
    4,2,
    7,0,
    5,1,
    3,3,
    3,5,
    1,5,
    5,3,
    4,0;
  
  const std::vector<std::vector<size_t>> memberships = {
    {1, 3, 8, 7},
    {2},
    {0, 6, 4, 5}
  };

  xtsp::ImplicitCompleteGraph<float> g(
    xyFull, std::make_shared<xtsp::Clustering>(xyFull.rows(), memberships));

  // the test subject
  auto gSuper = g.buildClusterMeans();
  ASSERT_EQ(gSuper.numVertices(), g.numClusters());
  const auto& xySuper = gSuper.getXy();
  EXPECT_FLOAT_EQ(xySuper(0,0), 4.5); EXPECT_FLOAT_EQ(xySuper(0,1), 1.5);
  EXPECT_FLOAT_EQ(xySuper(1,0), 7); EXPECT_FLOAT_EQ(xySuper(1,1), 0);
  EXPECT_FLOAT_EQ(xySuper(2,0), 2); EXPECT_FLOAT_EQ(xySuper(2,1), 4);
}