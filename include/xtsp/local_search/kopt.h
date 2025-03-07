#ifndef __XTSP_KOPT_H__
#define __XTSP_KOPT_H__

#include "xtsp/core/tour.h"
#include "xtsp/algorithm_utils/work_buffer.h"

namespace xtsp::algo
{
  template <typename CostTy>
  struct TwoOptQueryResults
  {
  public:
    CostTy improvement = 0; // this default value is important
    const size_t vA;                // index of vertex A
    size_t vC;                      // index of vertex C
  public:
    TwoOptQueryResults(size_t vA) : vA(vA)
    {
    }
    /// @brief Check if a new proposal make sense, if yes update 
    /// @param newImprovement the corresponding cost reduction
    /// @param vC_new
    /// @retval is the new improvement accepted?
    bool updateIfBetter(CostTy newImprovement, size_t vC_new)
    {
      if (newImprovement > this->improvement)
      {
        this->improvement = newImprovement;
        this->vC = vC_new;
        return true;
      }
      else
      {
        return false;
      }
    }
    bool isValid() const
    {
      return (improvement > 0);
    }
    // void apply(AbstractTour &tour, CostTy &originalCost, bool strict = false) const
    // {
    //   if (!isValid())
    //     throw std::invalid_argument("Trying to apply an invalid two-opt");
    //   tour.exchangeTwoEdges((size_t)vA, (size_t)vC, strict);
    // }
  };

  /// @brief For a given vertex A, find a 2-opt move
  ///        (either the first or best found.)
  ///
  /// Complexity: O(N) where N is the length of \p tour .
  /// To query if a (valid) move is found, simply query
  /// the method \p TwoOptQueryResults<CostTy>::isValid() .
  ///
  /// @param tour while not checked in the implementation, the tour
  ///        shall concern @p g and passes through @p vA exactly once)
  /// @param vA the ID of vertex A
  /// @param g the TSP instance
  /// @param firstImprovement false means "best improvement".
  /// @retval the two opt move, caller shall check if it's valid
  template <typename CostTy>
  TwoOptQueryResults<CostTy> find2OptMoveGivenA(
      const AbstractTour &tour, size_t vA, const AbstractCompGraph<CostTy> &g,
      bool firstImprovement);

  
  enum SweepMethod
  {
    kBitfieldTwoOptSweep, // this one is probably not as good as the priority queue approach.
    kPriorityTwoOptSweep,// turns out this is not new, similar to "Direct 2-Opt" in [GK10]
  };
  
  template <typename CostTy>
  struct TwoOptOutcome
  {
  public:
    CostTy improvement() const
    {
      return m_improvement;
    }
    size_t numMoves() const
    {
      return m_numMoves;
    }
    bool confirmedTwoOpt() const
    {
      return m_confirmedTwoOpt;
    }
    void update(CostTy extraImprovement, size_t extraMoves)
    {
      if (extraImprovement < 0 || m_confirmedTwoOpt)
        throw std::invalid_argument(
          "TwoOptOutcome<CostTy>::update expects an improvement to be +ve, "
          "and if it's already two-opt, you shouldn't call this again. "
          "Check if your code has bugs");
      m_improvement += extraImprovement;
      m_numMoves += extraMoves;
      m_confirmedTwoOpt = (extraMoves == 0);
    }
  protected: 
    CostTy m_improvement = 0;
    size_t m_numMoves = 0; // not just for book keeping 
    bool m_confirmedTwoOpt = false; 
  };

  ///
  /// @todo allow early termination?
  /// @todo try to avoid re-visit AB-CD and CD-AB twice
  /// @retval total cost improvement after one sweep of vertex A
  // already 2-opt???
  template <typename CostTy>
  class PriorityTwoOptFinder
  {
  public:
    // only memory allocation 
    PriorityTwoOptFinder(const AbstractTour &tour, const AbstractCompGraph<CostTy> &g);

    /// @param[out] numMovesDone after the call, this will be the amount of performed 2-opt moves
    TwoOptOutcome<CostTy> solve(
      AbstractTour &tour, const AbstractCompGraph<CostTy> &g,
      size_t maxNumSweeps = 10,
      bool firstImprovement = true);

    /// @param[out] numMovesDone (need not be 0), this function will increment it by the amount of moves
    TwoOptOutcome<CostTy> tryOneSweep2Opts(
      AbstractTour &tour, const AbstractCompGraph<CostTy> &g,
      bool firstImprovement = true);
  protected:
    void updateForNextSweep(
      const AbstractTour& tour, const AbstractCompGraph<CostTy> &g);
    
    // due to how C++ std::vector is typically implemented, 
    // it is more efficient to pop from back O(1) instead of the front O(N).
    // So we will arrange them in descending edge cost AB
    std::vector<std::pair<size_t, CostTy>> m_vAandCostAB;

    // a light-weight O(1) book-keeping to avoid
    // repeating AB-CD and CD-AB which is identical.
    // After all we do NOT maintain the edge costs values in 
    // m_vAandCostAB. So likely we will just be wasting resources
    // if we don't skip vertex C.
    std::vector<bool> m_skip; 
  };
}

#endif