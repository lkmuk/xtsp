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

  for (size_t i = 0; i < costMatrix.rows(); ++i)
    for (size_t j = 0; j < costMatrix.cols(); ++j)
      EXPECT_EQ(g.getEdgeCost(i,j), this->costMatrix(i,j)) 
        << "wrong edge cost value at (i,j) = (" << i << ", " << j << ")";
}

TEST_F(SimpleCompleteGraphInt, symmetricCostQuery)
{
  xtsp::CompleteGraph<int> g (true, this->costMatrix);

  EXPECT_TRUE(g.isSymmetric());
  EXPECT_FALSE(g.isClustered());

  for (size_t i = 0; i < costMatrix.rows(); ++i)
    for (size_t j = 0; j < costMatrix.cols(); ++j)
      EXPECT_EQ(g.getEdgeCost(i,j), g.getEdgeCost(j,i)) 
        << "not symmetric at (i,j) = (" << i << ", " << j << ")";

  // copied the lower-triangular?
  for (size_t i = 0; i < costMatrix.rows(); ++i)
    for (size_t j = 0; j < i; ++j)
      EXPECT_EQ(g.getEdgeCost(i,j), this->costMatrix(i,j)) 
        << "wrong edge cost value at (i,j) = (" << i << ", " << j << ")";
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
  xtsp::ImplicitCompleteGraph<float, 3> g (this->xyz);

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