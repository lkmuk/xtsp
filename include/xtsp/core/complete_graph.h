#ifndef __XTSP_CORE_COMPLETE_GRAPH_H__
#define __XTSP_CORE_COMPLETE_GRAPH_H__

#include "xtsp/core/clustering.h"

#include <memory>
#include <Eigen/Core>

// round a fp number to the nearest integer
// as suggested in TSPLIB
#define nint(x) ((int) (x+0.5))

namespace xtsp
{
  /// @brief Interfaces for a weighted complete graph.
  /// The graph can be either directed (for asymmetric TSP) or 
  /// undirected graph (for symmetric TSP).
  /// Optionally, you can declare graph partitioning (useful in GTSP).
  /// @tparam CostTy should be an unsigned integer/floating-point dtype
  template <typename CostTy>
  class AbstractCompGraph
  {
  public:
    virtual bool isSymmetric() const = 0;
    virtual CostTy getEdgeCost(size_t from, size_t to) const = 0;
    virtual size_t numVertices() const = 0;

    virtual size_t numClusters() const final;
    virtual bool isClustered() const final;
    virtual const std::shared_ptr<Clustering> getClusteringInfo() const final
    {
      return m_clustering;
    }

  protected:
    // how the vertices are clustered (nullptr means no clustering)
    const std::shared_ptr<Clustering> m_clustering = nullptr;
    
    // "protected" so that nobody can accidentally instantiate
    // this Abstract class
    AbstractCompGraph(const std::shared_ptr<Clustering> clustering = nullptr);

    // to remove Wdelete-non-abstract-non-virtual-dtor warning
    virtual ~AbstractCompGraph() = default;
  };

  /// @brief Weighted, complete graph with explicit representation of edge costs.
  /// @see AbstractCompGraph<CostTy>
  template <typename CostTy = int>
  class CompleteGraph : public AbstractCompGraph<CostTy>
  {
  public:
    /// @param isSymmetric 
    /// @param edgeCosts a NxN matrix. 
    ///         If it is symmetric, you can just fill out the lower triangle part.
    /// @param clustering 
    CompleteGraph(
      bool isSymmetric, 
      const Eigen::Matrix<CostTy, -1, -1>& edgeCosts, 
      const std::shared_ptr<Clustering> clustering = nullptr);

    virtual bool isSymmetric() const override;
    
    virtual CostTy getEdgeCost(size_t from, size_t to) const override
    {
      return m_mat(from, to);
    }
    virtual size_t numVertices() const override
    {
      return m_mat.rows();
    }

  protected:
    bool m_symmetric;
    // edge cost from "row Id" to "column Id"
    Eigen::Matrix<CostTy, -1, -1> m_mat; 
  };

  /// @brief Implicitly represent the N^2 edges where N is the number of vertices.
  /// For really large point set (say in millions),
  /// Most computers won't have enough RAM to store the whole cost matrix explicitly.
  /// So we create this class.
  ///
  /// We denote nDim as the number of dimension of each point.
  ///
  /// @tparam CostTy the base datatype for representing edge costs.
  ///         It should be a floating-point type!
  /// 
  ///
  /// @todo support explicitly pre-computation of the edge cost if N is small enough
  /// 
  template <typename CostTy = float> 
  class ImplicitCompleteGraph : public AbstractCompGraph<CostTy>
  {
  public:
    /// @param xy the list of points (N x nDim)
    /// @param clustering
    /// @param normType 2 -- Euclidean, 1 -- Manhattan, 0 -- maxNorm
    ImplicitCompleteGraph(
      const Eigen::Matrix<CostTy, -1, -1>& xy, 
      const std::shared_ptr<Clustering> clustering = nullptr,
      int normType = 2);

    /**
     * @brief Load a TSP or GTSP instance from a (G)TSPLIB-formatted file
     * 
     * The implementation is smart enough to deduce whether it's
     * a TSP and or GTSP instance.
     * 
     * @param fpath file path
     */
    static ImplicitCompleteGraph<CostTy> loadFromTsplibFile(
      const std::string& fpath);

    /// @brief the dimension of each data point (typically 2 or 3)
    size_t nDim() const;

    virtual bool isSymmetric() const override
    {
      return true;
    }
    virtual CostTy getEdgeCost(size_t from, size_t to) const override
    {
      switch (m_normType)
      {
        case 2:
          return (m_xy.row(from) - m_xy.row(to)).norm();
        case 1:
          return (m_xy.row(from) - m_xy.row(to)).array().abs().sum();
        default: // maxNorm, aka. L_infinty
          return (m_xy.row(from) - m_xy.row(to)).array().abs().maxCoeff();
      }
    }
    virtual size_t numVertices() const override
    {
      return m_xy.rows();
    }

    CompleteGraph<int> explicitize(float scale) const;

    /// @brief Construct super-graph where 
    /// each vertex corresponds to a cluster of this graph.
    /// Use case: local-global GTSP meta-heuristics.
    /// Here, simply use arithmetic averaging to determine 
    /// each cluster's (heuristic) "location".
    ///
    /// @throw if the graph is not clustered
    ImplicitCompleteGraph<CostTy> buildClusterMeans() const;

    const Eigen::Matrix<CostTy, -1, -1>& getXy() const;

  protected:
    Eigen::Matrix<CostTy, -1, -1> m_xy;
    int m_normType;
  };
} // namespace

#endif