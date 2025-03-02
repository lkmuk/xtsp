#include "xtsp/local_search/gtsp_only.h"

// Eliminate some logging at compile time
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_INFO
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>
#include <spdlog/fmt/bundled/ranges.h>
#include <cassert>

namespace xtsp::algo
{
  
  template <typename CostTy>
  GtspClusterOptimizer<CostTy>::GtspClusterOptimizer(size_t hint)
    : DynProgArena<CostTy>(hint)
  {
  }

  template <typename CostTy>
  CostTy GtspClusterOptimizer<CostTy>::solve(
    const xtsp::GeneralizedTour<CostTy>& tour, 
    const xtsp::AbstractCompGraph<CostTy>& graph, 
    size_t cutCluster,
    std::vector<size_t>& overallBestVertexSeq)
  {
    if (graph.getClusteringInfo() != tour.getClusteringInfo())
    {
      throw std::invalid_argument(
        "mismatch clustering info in the given graph and tour");
    }
    const std::shared_ptr<Clustering> clustering = graph.getClusteringInfo();
    const size_t numClusters = clustering->numClusters();
    if (cutCluster >= numClusters)
    {
      auto oldVal = cutCluster;
      cutCluster = cutCluster%numClusters;
      SPDLOG_WARN(
        "Reduced cutCluster fron {:d}, which was too large, to {:d}.",
        oldVal, cutCluster);
    }

    // only for debug build
    assert(tour.numClusters() == numClusters);
    assert(tour.numVertices() == graph.numVertices());

    // initialize the outputs 
    // (overall among the optimal solutions of all cut vertices)
    CostTy overallBestCost = std::numeric_limits<CostTy>::max();
    overallBestVertexSeq.clear();
    overallBestVertexSeq.resize(
      numClusters, std::numeric_limits<size_t>::max());

    SPDLOG_INFO("CO begins: cutCluster size = {:d})",
      clustering->getClusterSize(cutCluster));
    SPDLOG_DEBUG("cutCluster ID = {:d}, {:d} in total", cutCluster, numClusters);

    size_t rankCutCluster = tour.findClusterRankById(cutCluster);

    /// ideally this outer loop is iterated just once
    for (const size_t cutVertex : clustering->getMembers(cutCluster))
    {
      /// initialize
      this->clearBuf();
      // notice it's the number of vertices not the number of clusters M
      this->resizeBuf(graph.numVertices()); 

      /// backward-pass (start from the last row):
      ///       clusterSeq[rankCutCluster]
      ///   <-- clusterSeq[rankCutCluster + 1]
      ///   <-- clusterSeq[rankCutCluster + 2]
      ///   <-- ...
      ///   <-- clusterSeq[rankCutCluster + numClusters-1]
      ///
      ///@todo benchmark runtime-performance
      for (size_t i = 1; i < numClusters; ++i)
      {
        size_t rankPos = rankCutCluster + (numClusters - i);
        rankPos = rankPos%numClusters;
        size_t thisClusterId = tour.getClusterIdByRank(rankPos);
        const std::vector<size_t>& thisCluster = clustering->getMembers(thisClusterId);
        for (const size_t vertex : thisCluster)
        {
          // a terminal state (after that, the tour goes to cutVertex again)
          if (i==1) 
          {
            SPDLOG_DEBUG("terminal vertex: n = {:d}, m = {:d}, rank: {:d}", 
              vertex, thisClusterId, rankPos);
            this->bestNextVertex[vertex] = cutVertex;
            // the terminal cost of this state
            this->costToGo[vertex] = graph.getEdgeCost(vertex, cutVertex);
          }
          else
          {
            SPDLOG_DEBUG("intermediate vertex: n = {:d}, m = {:d}, rank: {:d}", 
              vertex, thisClusterId, rankPos);
            size_t nextClusterId = tour.getClusterIdByRank((rankPos+1)%numClusters);
            const std::vector<size_t>& nextCluster = clustering->getMembers(nextClusterId); 
            this->backpassStep(
              vertex, nextCluster, 
              [& graph](size_t from, size_t to){return graph.getEdgeCost(from,to);});
          }
        }
      }

      /// the last backward pass (the step from cutVertex to the next cluster's)
      size_t nextClusterId = tour.getClusterIdByRank((rankCutCluster+1)%numClusters);
      const std::vector<size_t>& nextCluster = clustering->getMembers(nextClusterId); 
      this->backpassStep(
        cutVertex, nextCluster,
        [& graph](size_t from, size_t to){return graph.getEdgeCost(from,to);});

      /// The forward-pass.
      /// We only need to do it if this cutVertex leads to a new best.
      ///
      /// Alternatively, we could also just copy the costToGo and bestNextVertex
      /// to somewhere
      /// and defer the forward pass until the end.
      if (this->costToGo[cutVertex] < overallBestCost)
      {
        overallBestCost = this->costToGo[cutVertex];
        SPDLOG_INFO(
            "new best: tour cost = {}, cutVertex = {:d}", 
            overallBestCost, cutVertex);
        SPDLOG_DEBUG(
            "Cost-to-go table:\n{}\n"
            "Best next move table:\n{}\n",
            fmt::join(this->costToGo, " "),
            fmt::join(this->bestNextVertex, " "));
        // This will automatically follow the same cluster sequence.
        // Here, overallBestVertexSeq already has the right size.
        overallBestVertexSeq[rankCutCluster] = cutVertex;
        for (size_t delta = 1, thisVertex_ = cutVertex; delta < numClusters; ++delta)
        {
          size_t nextBestVertex = this->bestNextVertex[thisVertex_];
          overallBestVertexSeq[(rankCutCluster+delta)%numClusters] = nextBestVertex;
          thisVertex_ = nextBestVertex;
        }
      }
    }
    return overallBestCost;
  }

  template <typename CostTy>
  void GtspClusterOptimizer<CostTy>::improve(
    xtsp::GeneralizedTour<CostTy> &tour, 
    const xtsp::AbstractCompGraph<CostTy> &graph, 
    size_t cutCluster)
  {
    // during the solve, we don't need to use tour.m_seq.
    // So we can safely "reuse" it as the output buffer for
    // the `solve` call. 
    // This guarantees that there is no memory allocation
    // during this `improve` routine.
    std::vector<size_t> &newTour = tour.getSeqMutableRef__();
    CostTy newCost = this->solve(tour, graph, cutCluster, newTour);
    tour.getCostMutableRef() = newCost;
    
    // // In the case of Cluster Optimization, 
    // // no need to update the cache
    // tour.updateCachedClusterSeq();
  }


  // explicit template instantiation
  template class GtspClusterOptimizer<float>;
  template class GtspClusterOptimizer<int>;

} // namespace xtsp::algo