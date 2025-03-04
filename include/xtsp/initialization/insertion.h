#ifndef __XTSP_TOUR_CONSTRUCTION_INSERTION_H__
#define __XTSP_TOUR_CONSTRUCTION_INSERTION_H__

#include "xtsp/core/tour.h"

namespace xtsp::algo
{
  template <typename CostTy>
  Tour<CostTy> farthestInsertion(const AbstractCompGraph<CostTy>& g);
}

#endif