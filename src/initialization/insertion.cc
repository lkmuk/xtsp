#include "xtsp/initialization/insertion.h"

#include <set>

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_WARN
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/ranges.h>

namespace xtsp::algo
{
  template <typename CostTy>
  Tour<CostTy> farthestInsertion(const AbstractCompGraph<CostTy>& g)
  {
    // a straight-forward implementation 

    // initialize a partial Hamiltonian tour of length 2
    const size_t numV = g.numVertices(); // aka N
    assert(numV >= 2);
    std::vector<size_t> pTour; // partial Hamiltonian tour: rank --> vertex ID
    pTour.emplace_back(0);
    pTour.emplace_back(numV-1);

    /// @todo should we use bitfield instead?
    std::set<size_t> verticesToAdd;
    for (size_t i = 0; i < numV; ++i)
    {
      // if not in pTour
      if (std::find(pTour.cbegin(), pTour.cend(), i) == pTour.cend())
        verticesToAdd.emplace_hint(verticesToAdd.end(), i);
    }
    pTour.reserve(numV);
    SPDLOG_DEBUG("Initialized a partial tour: {}", fmt::join(pTour, " "));


    // iterate for N-1 times
    size_t vPicked;
    while (!verticesToAdd.empty())
    {
      // which one to insert? (farthest from the pTour)
      CostTy maxDistFromTour = std::numeric_limits<CostTy>::min();
      for (const auto vX : verticesToAdd)
      {
        // compute the distance of X from the partial tour
        CostTy distTour2X = std::numeric_limits<CostTy>::max();
        for (const auto vExisting : pTour)
        {
          CostTy d = g.getEdgeCost(vExisting,vX);
          if (d < distTour2X)
          {
            distTour2X = d;
          }
        }

        // update vPicked and maxDistFromTour
        if (distTour2X > maxDistFromTour)
        {
          maxDistFromTour = distTour2X;
          vPicked = vX;
        }
      }
    
      // where to insert? (greedy: minimize insertion cost)
      CostTy minInsertionCost = std::numeric_limits<CostTy>::max();
      size_t whereToInsert; // note: vPicked will be inserted before whereToInsert
      for (size_t i = 0; i < pTour.size(); ++i) 
      {
        size_t vWhere = pTour[i];
        size_t vPrev = (i == 0) ? pTour.back() : pTour[i-1];
        CostTy insertionCost = g.getEdgeCost(vPrev, vPicked) + g.getEdgeCost(vPicked, vWhere) - g.getEdgeCost(vPrev, vWhere);
        if (insertionCost < minInsertionCost)
        {
          minInsertionCost = insertionCost;
          whereToInsert = i; // rank instead of vertex ID
        }
      }
      
      // do the insertion
      SPDLOG_DEBUG("Inserting vertex {:d} at rank {:d}", vPicked, whereToInsert);
      pTour.insert(pTour.begin()+whereToInsert, vPicked);
      verticesToAdd.erase(vPicked);
    } 

    Tour<CostTy> tour(pTour, numV);
    tour.evalTour(g);
    return tour;
  }

  template Tour<float> farthestInsertion(const AbstractCompGraph<float>&);
  template Tour<int> farthestInsertion(const AbstractCompGraph<int>&);
}