#ifndef __XTSP_LOCAL_SEARCH_DP_H__
#define __XTSP_LOCAL_SEARCH_DP_H__

#include <cstddef>
#include <vector>
#include <functional>

namespace xtsp::algo
{
  template <typename CostTy>
  struct DynProgArena
  {
    std::vector<size_t> bestNextVertex;
    std::vector<CostTy> costToGo;

    /// @brief pre-allocate Dynamic Programming work memory
    ///        instead of global variables!
    /// @param[in] hintNumVertices just a hint, the work data can expand during runtime.
    DynProgArena(size_t hintNumVertices)
    {
      bestNextVertex.reserve(hintNumVertices);
      costToGo.reserve(hintNumVertices);
    }

    /// @brief just reset while keeping the allocated memory
    void clearBuf()
    {
      bestNextVertex.clear();
      costToGo.clear();
    }

    /// @brief resize and zero initialize the data
    void resizeBuf(size_t newSize)
    {
      bestNextVertex.resize(newSize, 0);
      costToGo.resize(newSize, 0);
    }

    /// @brief one-step look-ahead update
    /// @pre all possible next vertices's costToGo values are already populated.
    void backpassStep(
      size_t fromVertex,
      const std::vector<size_t>& allPossibleNextVertices,
      std::function<CostTy(size_t,size_t)> edgeCostFnc)
    {
      size_t bestQNextVertex = std::numeric_limits<size_t>::max();
      CostTy bestQVal = std::numeric_limits<CostTy>::max();
      for (const size_t nextVertex : allPossibleNextVertices)
      {
        CostTy qVal = edgeCostFnc(fromVertex, nextVertex) 
                      + this->costToGo[nextVertex];
        if (qVal < bestQVal)
        {
          bestQNextVertex = nextVertex;
          bestQVal = qVal;
        }
      }
      this->costToGo[fromVertex] = bestQVal;
      this->bestNextVertex[fromVertex] = bestQNextVertex;
    }
  };
}

#endif