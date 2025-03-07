#include "xtsp/local_search/kopt.h"
#include "xtsp/core/tour_alternatives.h"
#include "xtsp/local_search/gtsp_only.h"

#include <gtest/gtest.h>
#include <iostream>
#include <spdlog/spdlog.h>

class TwoOptHamTour : public testing::Test
{
protected:
  Eigen::MatrixX2f m_pts{
    // {0,0},
    // {7,0},
    // {0,2},
    // {3,3},
    // {7,3},
    // {0,6},
    // {3,6},
    // {3, 10},
    // {5, 10}   // optimal dataset for this test case 0, 1, 4, 3, 5, 8, 7, 6, 2
    // a better test case (e.g., 3-5 crosses both 7-2 and 2-6)
    {0,0},
    {4,0},
    {1,2},
    {3,2},
    {4,2},
    {0,4},
    {2,4},
    {4,4},
    {2,6}
  };
  xtsp::ImplicitCompleteGraph<float> m_graph {m_pts};
  xtsp::AdjTabTour m_tour {std::vector<size_t>({
    0, 1, 4, 8, 7, 2, 6, 3, 5
  })};
  xtsp::AdjTabTour m_tourOpt {std::vector<size_t>({
    // 0, 1, 4, 7, 3, 8, 6, 2, 5 // wrong!
    0, 1, 3, 4, 7, 8, 6, 5, 2
  })};
  void SetUp() override
  {
    // check if the test has been implemented correctly
    ASSERT_EQ(m_graph.numVertices(), 9);
    ASSERT_TRUE(m_tour.isHamiltonian());
  }
};

TEST_F(TwoOptHamTour, cannotFind2OptGivenA)
{
  spdlog::set_level(spdlog::level::debug);

  size_t vA = 0;
  spdlog::debug("First improvement");
  auto resFirstImprov = xtsp::algo::find2OptMoveGivenA(m_tour, vA, m_graph, true);
  EXPECT_FALSE(resFirstImprov.isValid());

  spdlog::debug("Best improvement");
  auto resBestImprov = xtsp::algo::find2OptMoveGivenA(m_tour, vA, m_graph, false);
  EXPECT_FALSE(resBestImprov.isValid());
}

TEST_F(TwoOptHamTour, canFind2OptGivenA)
{
  spdlog::set_level(spdlog::level::debug);
  size_t vA = 2;
  size_t expectedVertexC = 3; // based on our test case

  spdlog::debug("First improvement");
  auto resFirstImprov = xtsp::algo::find2OptMoveGivenA(m_tour, vA, m_graph, true);
  EXPECT_TRUE(resFirstImprov.isValid());
  EXPECT_EQ(resFirstImprov.vC, expectedVertexC);

  spdlog::debug("Best improvement");
  auto resBestImprov = xtsp::algo::find2OptMoveGivenA(m_tour, vA, m_graph, false);
  EXPECT_TRUE(resBestImprov.isValid());
  EXPECT_EQ(resFirstImprov.vC, expectedVertexC);
}

TEST_F(TwoOptHamTour, solveFirstImprovMax10Sweeps)
{
  spdlog::set_level(spdlog::level::debug);
  auto oldTourCost = xtsp::evalTour(m_tour, m_graph);
  auto minTourCost = xtsp::evalTour(m_tourOpt, m_graph);

  xtsp::algo::PriorityTwoOptFinder<float> solver(m_tour, m_graph);
  auto res = solver.solve(m_tour, m_graph, 10, true);
  EXPECT_TRUE(res.confirmedTwoOpt());
  EXPECT_TRUE(res.improvement() > 0);
  EXPECT_FLOAT_EQ(oldTourCost - res.improvement(), minTourCost);
  m_tour.print(std::cout);
  std::cout << std::endl;
}

TEST_F(TwoOptHamTour, solveBestImprovMax10Sweeps)
{
  spdlog::set_level(spdlog::level::debug);
  auto oldTourCost = xtsp::evalTour(m_tour, m_graph);
  auto minTourCost = xtsp::evalTour(m_tourOpt, m_graph);

  xtsp::algo::PriorityTwoOptFinder<float> solver(m_tour, m_graph);
  auto res = solver.solve(m_tour, m_graph, 10, false);
  EXPECT_TRUE(res.confirmedTwoOpt());
  EXPECT_TRUE(res.improvement() > 0);
  EXPECT_FLOAT_EQ(oldTourCost - res.improvement(), minTourCost);
  m_tour.print(std::cout);
  std::cout << std::endl;
}
/// @todo test partial Hamiltonian tour
