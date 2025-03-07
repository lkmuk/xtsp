// #define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#include <spdlog/spdlog.h>

#include "xtsp/local_search/kopt.h"

namespace xtsp::algo
{
  template <typename CostTy>
  TwoOptQueryResults<CostTy> find2OptMoveGivenA(
      const AbstractTour &tour, size_t vA, const AbstractCompGraph<CostTy> &g,
      bool firstImprovement)
  {
    if (tour.size() < 4)
      throw std::invalid_argument("Your tour is too short");
    if (!g.isSymmetric())
      throw std::invalid_argument(
          "Currently the 2-opt implementation doesn't support assymmetric TSP yet");
    TwoOptQueryResults<CostTy> result(vA);
    const auto vB = tour.next(vA);
    const CostTy cAB = g.getEdgeCost(vA, vB);
    auto vC = tour.next(vB);
    // notice the termination condition 
    // (which guarantee each iteration ABCD is valid);
    for (size_t ii = 0; ii < tour.size() - 3; ++ii)
    {
      auto vD = tour.next(vC);
      // our for-loop should never let vD == vA
      SPDLOG_DEBUG(
        "ii = {:d}: A,B,C,D = {:d},{:d},{:d},{:d} (tour: {})",
        ii, vA, vB, vC, vD, tour.print());
      assert(vD != vA);

      CostTy oldComponent = cAB + g.getEdgeCost(vC, vD);
      CostTy newComponent = g.getEdgeCost(vA, vC) + g.getEdgeCost(vB, vD);
      // we assume flipping either segment BC or AD
      // has no impact on their respective tour cost component
      // (we can take care of it later if we want to support ATSP)
      CostTy improvement = oldComponent - newComponent;
      bool isAccepted = result.updateIfBetter(improvement, vC);
      if (isAccepted)
      {
        SPDLOG_DEBUG("accepted the move with improvement: {}", improvement);
      }
      if (isAccepted && firstImprovement)
      {
        return result;
      }

      // update for the next iter
      vC = tour.next(vC);
    }
    // if not found, return a dummy;
    // otherwise, this corresponds to the best 2-opt move for this given vertex A
    return result;
  }

  template <typename CostTy>
  PriorityTwoOptFinder<CostTy>::PriorityTwoOptFinder(
      const AbstractTour &tour, const AbstractCompGraph<CostTy> &g)
  {
    if (tour.maxSize() != g.numVertices())
      throw std::invalid_argument(
          "PriorityTwoOptFinder ctor: tour is inconsistent with the graph");
    updateForNextSweep(tour, g);
  }

  template <typename CostTy>
  void PriorityTwoOptFinder<CostTy>::updateForNextSweep(
      const AbstractTour &tour, const AbstractCompGraph<CostTy> &g)
  {
    // we allow the tour to be partially-Hamiltonian (e.g., in generalized TSP )
    m_skip = std::vector<bool>(tour.maxSize(), false);
    m_vAandCostAB.clear();
    m_vAandCostAB.reserve(tour.size());
    size_t vA = tour.getDepotId();
    for (size_t rank = 0; rank < tour.size(); ++rank)
    {
      auto vB = tour.next(vA);
      m_vAandCostAB.emplace_back(
        std::make_pair<size_t &, CostTy>(vA, g.getEdgeCost(vA, vB)));
      vA = vB;
    }
    // sorting (in descending!!! edge cost AB)
    std::sort(m_vAandCostAB.begin(), m_vAandCostAB.end(),
              [](const std::pair<size_t, CostTy> &lhs, std::pair<size_t, CostTy> &rhs)
              {
                return lhs.second > rhs.second;
              });
  }

  template <typename CostTy>
  TwoOptOutcome<CostTy> PriorityTwoOptFinder<CostTy>::tryOneSweep2Opts(
      AbstractTour &tour, const AbstractCompGraph<CostTy> &g,
      bool firstImprovement)
  {
    TwoOptOutcome<CostTy> outcome;
    updateForNextSweep(tour, g);
    for (const auto &[vA, preComputedCostAB] : m_vAandCostAB)
    {
      // the preComputedCostAB is more for prioritizing vA,
      // during the actual twoOpt search, we will lead to we evaluate
      // them because the precomputed cost value can be outdated
      if (m_skip[vA])
        continue;
      auto res = find2OptMoveGivenA<CostTy>(tour, vA, g, firstImprovement);
      assert(res.vA == vA);
      m_skip[vA] = true;
      if (res.isValid())
      {
        /// @todo confirm empirically if this speeds up the improvement
        m_skip[res.vC] = true;

        outcome.update(res.improvement, 1);
        /// in this Priority-based method, it might be better to 
        /// specify the sequence. @todo evidence
        SPDLOG_DEBUG("Perform a two-opt move: A = {:d}, C = {:d}", res.vA, res.vC);
        SPDLOG_DEBUG("Tour (currently): {}", tour.print());
        tour.exchangeTwoEdges(res.vA, res.vC, true);
        SPDLOG_DEBUG("Tour (new)      : {}", tour.print());
      }
    }
    
    return outcome;
  }

  template <typename CostTy>
  TwoOptOutcome<CostTy> PriorityTwoOptFinder<CostTy>::solve(
      AbstractTour &tour, const AbstractCompGraph<CostTy> &g,
      size_t maxNumSweeps,
      bool firstImprovement)
  {
    {
      std::string displayName;
      if (firstImprovement)
        displayName = "first";
      else
        displayName = "best";
      SPDLOG_INFO(
        "priority 2-opt: {}-improvement, max. {:d} sweep(s)",
        displayName, maxNumSweeps);
    }
    if (maxNumSweeps == 0)
      SPDLOG_WARN("priority 2-opt: Ignoring no-op request");
    TwoOptOutcome<CostTy> overallResult; // overall across all sweeps so far
    // CostTy improvementLastSweep = 1e3; // a dummy value
    for (size_t totNumSweeps = 0; totNumSweeps < maxNumSweeps; ++totNumSweeps)
    {
      SPDLOG_INFO("sweep {:d} begins", totNumSweeps+1);

      TwoOptOutcome<CostTy> sweepRes = tryOneSweep2Opts(
        tour, g, firstImprovement);
      
      SPDLOG_INFO(
          "sweep {:d} : further improved by {} using {:d} moves", 
          totNumSweeps+1, sweepRes.improvement(), sweepRes.numMoves());
      overallResult.update(sweepRes.improvement(), sweepRes.numMoves());
      /// @todo early termination if the the improvement is 
      /// very small for several sweeps in a row
      if (sweepRes.numMoves() != 0) // i.e., > 0
      {
        assert(sweepRes.improvement() > 0);
      }
      else // i.e. (numMovesThisSweep == 0)
      {
        SPDLOG_INFO(
          "no move found, so 2-opt is confirmed", 
          totNumSweeps+1);
        return overallResult;
      }

      // improvementLastSweep = sweepRes.improvement;
    }
    return overallResult;
  }

  /**********************************
    Explict template instantiation
   **********************************/
  template TwoOptQueryResults<int> find2OptMoveGivenA(
      const AbstractTour &tour, size_t vA, const AbstractCompGraph<int> &g,
      bool firstImprovement);
  template TwoOptQueryResults<float> find2OptMoveGivenA(
      const AbstractTour &tour, size_t vA, const AbstractCompGraph<float> &g,
      bool firstImprovement);
  template class PriorityTwoOptFinder<float>;
  template class PriorityTwoOptFinder<int>;

}