#include "xtsp/initialization/insertion.h"

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_WARN
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/ranges.h>

#include "xtsp/algorithm_utils/work_buffer.h"

namespace xtsp::algo
{
  template <typename CostTy>
  PermTour farthestInsertion(const AbstractCompGraph<CostTy> &g, size_t vFirstPick)
  {
    // a straight-forward implementation

    if (vFirstPick >= g.numVertices())
    {
      std::string errMsg = fmt::format(
          "Invalid first pick ID {:d} (should be 0 <= pick < N = {:d})",
          vFirstPick, g.numVertices());
      SPDLOG_ERROR(errMsg);
      throw std::invalid_argument(errMsg);
    }

    // initialize a partial Hamiltonian tour of length 2
    const size_t numV = g.numVertices(); // aka N
    assert(numV >= 2);
    std::vector<size_t> pTour; // partial Hamiltonian tour: rank --> vertex ID
    pTour.emplace_back(vFirstPick);

    xtsp::internal::WorkBuffer verticesToAdd(numV);
    for (auto i : pTour)
    {
      verticesToAdd.markDone(i);
    }
    pTour.reserve(numV);
    SPDLOG_DEBUG("Initialized a partial tour: {}", fmt::join(pTour, " "));

    // iterate for N-1 times
    size_t vPicked;
    while (!verticesToAdd.isEmpty())
    {
      // which one to insert? (farthest from the pTour)
      CostTy maxDistFromTour = std::numeric_limits<CostTy>::min();
      /// @todo if we are obessessed with runtime,
      ///     consider looping over the min and max of verticesToAdd instead
      for (size_t vX = verticesToAdd.estimateTodoIdMin();
           vX <= verticesToAdd.estimateTodoIdMax(); ++vX)
      {
        if (!verticesToAdd.isDone(vX))
          continue;
        // compute the distance of X from the partial tour
        CostTy distTour2X = std::numeric_limits<CostTy>::max();
        for (const auto vExisting : pTour)
        {
          CostTy d = g.getEdgeCost(vExisting, vX);
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
      size_t whereToInsert; // note: vPicked will be inserted before whereToInsert
      if (pTour.size() == 1)
      {
        whereToInsert = 1; // or 0
      }
      else
      {
        CostTy minInsertionCost = std::numeric_limits<CostTy>::max();
        for (size_t i = 0; i < pTour.size(); ++i)
        {
          size_t vWhere = pTour[i];
          size_t vPrev = (i == 0) ? pTour.back() : pTour[i - 1];
          CostTy insertionCost = g.getEdgeCost(vPrev, vPicked) + g.getEdgeCost(vPicked, vWhere) - g.getEdgeCost(vPrev, vWhere);
          if (insertionCost < minInsertionCost)
          {
            minInsertionCost = insertionCost;
            whereToInsert = i; // rank instead of vertex ID
          }
        }
      }

      // do the insertion
      SPDLOG_DEBUG("Inserting vertex {:d} at rank {:d}", vPicked, whereToInsert);
      pTour.insert(pTour.begin() + whereToInsert, vPicked);
      verticesToAdd.markDone(vPicked);
    }

    PermTour tour(pTour, numV);
    // tour.evalTour(g);
    return tour;
  }

  template PermTour farthestInsertion<float>(const AbstractCompGraph<float> &, size_t);
  template PermTour farthestInsertion<int>(const AbstractCompGraph<int> &, size_t);
}