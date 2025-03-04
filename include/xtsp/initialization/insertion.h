#ifndef __XTSP_TOUR_CONSTRUCTION_INSERTION_H__
#define __XTSP_TOUR_CONSTRUCTION_INSERTION_H__

#include "xtsp/core/tour.h"

namespace xtsp::algo
{
  /// @brief construct a valid Hamiltonian tour via farthest insertion
  /// @param vFirstPick ID of the first vertex to be added to the empty tour 
  template <typename CostTy>
  Tour<CostTy> farthestInsertion(const AbstractCompGraph<CostTy>& g, size_t vFirstPick = 0);
}

#endif