#include "xtsp/local_search/gtsp_only.h"

#include <gtest/gtest.h>
#include <iostream>
#include <spdlog/spdlog.h>

class ClusterOptimizerToyExample : public testing::Test
{
protected:
  static void add3x3SquareGrid(Eigen::MatrixX2f& out, const Eigen::RowVector2f& center)
  {
    auto head = out.rows();
    out.conservativeResize(head+9, 2);
    for (int col = -1; col <= 1; ++col)
      for (int row = -1; row <= 1; ++row)
      {
        out.row(head++) = Eigen::RowVector2f({col, row})+center;
      }
  }
  virtual void SetUp () override
  {
    size_t trueNumVertices = 36;
    
    // construct the points
    Eigen::Matrix<float, -1, 2> xy(trueNumVertices, 2);
    xy.resize(0, 2); // after reserving the memory, shrink it again
    add3x3SquareGrid(xy, {0,0}); // cluster 0
    add3x3SquareGrid(xy, {5,5}); // cluster 1
    add3x3SquareGrid(xy, {0,5}); // cluster 2

    xy.conservativeResize(trueNumVertices, 2);
    // std::cout << "Point set (clusters 3, 4, 5 remain to be initialized)\n";
    // std::cout << xy << std::endl;

    // cluster 3
    xy.row(27) = Eigen::RowVector2f({3, 6});
    xy.row(28) = xy.row(27) + Eigen::RowVector2f({1, -1});
    xy.row(29) = xy.row(27) + Eigen::RowVector2f({0, -2});
    xy.row(30) = xy.row(27) + Eigen::RowVector2f({0, -3});
    xy.row(31) = xy.row(27) + Eigen::RowVector2f({1, -3});

    // cluster 4
    xy.row(32) = Eigen::RowVector2f({4, -4});
    xy.row(33) = Eigen::RowVector2f({6, -3});
    xy.row(34) = Eigen::RowVector2f({4, -3});

    // cluster 5
    xy.row(35) = Eigen::RowVector2f({4, 1});

    // std::cout << "Point set (final)\n";
    // std::cout << xy << std::endl;

    // construct the cluster.
    // here a point = a GTSP city = a vertex
    std::vector<size_t> numPointsPerCluster {9, 9, 9, 5, 3, 1};
    m_clustering = std::make_shared<xtsp::Clustering>(
      xtsp::Clustering::cumsum(numPointsPerCluster));
    ASSERT_EQ(m_clustering->numVertices(), trueNumVertices) 
      << "Bug in the implementation of the test";
    
    // construct the graph
    m_graph = std::make_unique<xtsp::ImplicitCompleteGraph<float>>(
      xy, m_clustering
    );

    // construct a tour according to the said cluster sequence
    std::vector <size_t> initTour;
    // tour rank "tr"
    for (size_t tr = 0; tr < m_clustering->numClusters(); ++tr)
    {
      // all vertices among the cluster at this tour position
      auto members = m_clustering->getMembers(m_clusterSeq[tr]);
      // pick the vertex of the lowest index.
      initTour.emplace_back(members[0]);
    }
    m_tour = std::make_unique<xtsp::GeneralizedTour<float>>(
      initTour, m_clustering);
    for (size_t tr = 0; tr < m_clustering->numClusters(); ++tr)
    {
      ASSERT_EQ(m_clusterSeq[tr], m_tour->getClusterIdByRank(tr)) 
        << "Wrong test implementation: at rank = " 
        << tr << ", expect " << m_clusterSeq[tr] 
        << " but got " << m_tour->getClusterIdByRank(tr);
    }

    // float oldTourCost = m_tour->evalTour(*m_graph);
    // std::cout << "Original cost = " << oldTourCost << std::endl;
    xtsp::GeneralizedTour<float> optimizedTour(m_trueOptGenTour, m_clustering);
    float reevaluatedMinTourCost = optimizedTour.evalTour(*m_graph);
    EXPECT_FLOAT_EQ(reevaluatedMinTourCost, m_trueOptCost) 
      << "Buggy test implementation: Wrong computed optimal tour cost";
  }
  std::shared_ptr<xtsp::Clustering> m_clustering = nullptr;
  // we will somehow initialize it so that it satisfy the cluster sequence 
  std::unique_ptr<xtsp::GeneralizedTour<float>> m_tour = nullptr;
  std::unique_ptr<xtsp::ImplicitCompleteGraph<float>> m_graph = nullptr;
  const std::vector<size_t> m_clusterSeq {0, 4, 5, 1, 3, 2};
  // The expected generalized tour, in which the cluster sequence is kept
  const std::vector<size_t> m_trueOptGenTour {8, 34, 35, 9, 29, 24};
  const float m_trueOptCost = 5 + 4 + 3 + 2 + 1 + 3;
};

TEST_F(ClusterOptimizerToyExample, solveApiWorks)
{
  spdlog::set_level(spdlog::level::warn);
  auto cutCluster = m_clustering->evalWhichHasTheLeastVertices();

  std::vector<size_t> computedOptVertices;
  // System under test
  xtsp::algo::GtspClusterOptimizer<float> solver(m_clustering->numVertices());
  // local optimal cost (while fixing the cluster sequence)
  float computedOptCost = solver.solve(*m_tour, *m_graph, cutCluster, computedOptVertices);
  ASSERT_EQ(computedOptVertices.size(), m_clustering->numClusters());
  EXPECT_FLOAT_EQ(computedOptCost, m_trueOptCost);
  for (size_t rank = 0; rank < m_clustering->numClusters(); ++rank)
  {
    EXPECT_EQ(computedOptVertices[rank], m_trueOptGenTour[rank]) 
      << "at tour rank = " << rank;
  }
}

TEST_F(ClusterOptimizerToyExample, recommendedApi)
{
  spdlog::set_level(spdlog::level::warn);
  auto cutCluster = m_clustering->evalWhichHasTheLeastVertices();

  // System under test
  xtsp::algo::GtspClusterOptimizer<float> solver(m_clustering->numVertices());
  solver.improve(*m_tour, *m_graph, cutCluster);

  ASSERT_EQ(m_tour->size(), m_clustering->numClusters());
  EXPECT_FLOAT_EQ(m_tour->getCost(), m_trueOptCost);
  for (size_t rank = 0; rank < m_clustering->numClusters(); ++rank)
  {
    EXPECT_EQ(m_tour->getVertex(rank), m_trueOptGenTour[rank]) 
      << "at tour rank = " << rank;
  }
}

TEST_F(ClusterOptimizerToyExample, tryAllClustersAsCut)
{
  spdlog::set_level(spdlog::level::info);
  // Preallocate the work and output buffers
  std::vector<size_t> computedOptVertices;
  xtsp::algo::GtspClusterOptimizer<float> solver(m_clustering->numVertices());
  
  for (size_t cutCluster = 0; cutCluster < m_clustering->numClusters(); ++cutCluster)
  {
    // local optimal cost (while fixing the cluster sequence)
    float computedOptCost = solver.solve(*m_tour, *m_graph, cutCluster, computedOptVertices);
    ASSERT_EQ(computedOptVertices.size(), m_clustering->numClusters());
    EXPECT_FLOAT_EQ(computedOptCost, m_trueOptCost);
    for (size_t rank = 0; rank < m_clustering->numClusters(); ++rank)
    {
      EXPECT_EQ(computedOptVertices[rank], m_trueOptGenTour[rank]) 
        << "at tour rank = " << rank;
    }
  }
}

/// @todo test scalability