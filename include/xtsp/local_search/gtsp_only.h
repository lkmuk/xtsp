#ifndef __XTSP_GTSP_ONLY_ALGO_H__
#define __XTSP_GTSP_ONLY_ALGO_H__

#include "xtsp/core/complete_graph.h"
#include "xtsp/core/tour.h"
#include "xtsp/dyn_prog.h"

namespace xtsp
{
  namespace algo
  {
    /** 
     * @brief the so-called Cluster Optimization as per [GK10]
     * 
     * ## Quick overview
     * 
     * This local search technique was proposed as early as 
     * [FST97] (see Figure 5 of that paper.)
     * Given a (fixed) cluster visit sequence, it efficiently 
     * finds for each cluster the optimal city to visit for each cluster.
     *  
     *  * It solves the local search exactly.
     *  * Applicable to both symmetric and asymmetric GTSP
     *  * Used in [FST97] 
     *  * Used in Memetic Algorithm of [GK10].
     *  * can also be used for the **local-global heuristics** 
     *    for generalized TSP. See ???
     * 
     * @ref [FST97] Fischetti, Salazar-Gonzalez, & Toth. (1997). 
     *      A branch-and-cut algorithm for the symmetric generalized 
     *      traveling salesman problem. Operations Research, 45(3), 378-394.
     * @todo update more literature reference
     */
    template <typename CostTy>
    class GtspClusterOptimizer : public DynProgArena<CostTy>
    {
    public:
      /// @brief Preallocate the work buffer based on the hint
      /// @param reserveSize a hint on the maximum total number of vertices
      GtspClusterOptimizer(size_t reserveSize);

      /**
       * @brief solve a new Cluster Optimization problem
       * 
       * ## About the algorithm's efficiency
       * 
       * Let $M$ be the number of clusters,
       *     $N$ be the total number of vertices,
       *     $x$ be the number of vertices in \p cutCluster , and
       *     $s$ be the maximum number of vertices in a cluster.
       * 
       * The optimization is exact. Nonetheless it can be 
       * efficiently solved as $x$ Dynamic Programs (DP).
       * Each DP involves $M-1$ time-steps, which is a
       * well-known kind of DP, especially in Optimal Control:
       * it begins with a **backward pass/recursion**, 
       * followed by a forward recursion.
       * The most expensive part is the backward pass:
       * combining all $x$ DP (implemented sequentially)
       * 
       * * the space complexity is $O(N)$.
       * * the time complexity can also be expressed analytically.
       *   To simplify the expression, we state only the 
       *   worst-case complexity of $O(x * (M-1) * s^2)$,
       *   i.e., polynomial w.r.t. $s$ and linear w.r.t. $x$ and $M$.
       *   Compare it to exhaustive serach's $O(s^M)$.
       * 
       * ## An implementation detail
       * 
       * Conceptually, we create another graph with the vertices
       * in the cut cluster being "duplicated";
       * In the actual implementation, we don't 
       * copy the cut cluster. 
       * 
       * @param[in] tour 
       * @param[in] graph 
       * @param[in] cutCluster 
       *            The index of the cluster where we cut the generalized tour.
       *            For efficient execution you should specify a singleton cluster,
       *            i.e., with only one city. 
       *            You should also ensure it's within [0, ..., N-1].
       * @param[out] optimalTour The optimized generalized tour,
       *            which follows the same original cluster sequence 
       *             suggested in \p tour .
       * @retval the new tour cost after the cluster optimization
       * @sa \p xtsp::Clustering::evalWhichHasTheLeastVertices
       */
      CostTy solve(
        const xtsp::GeneralizedTour<CostTy> &tour, 
        const xtsp::AbstractCompGraph<CostTy> &graph, 
        size_t cutCluster,
        std::vector<size_t>& optimalTour);
    };

  } // namespace algo  
} // namespace xtsp

#endif
