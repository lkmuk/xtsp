#ifndef __XTSP_CORE_TOUR_H__
#define __XTSP_CORE_TOUR_H__

#include "xtsp/core/complete_graph.h"

#include <string_view>
#include <limits>
#include <type_traits>
#include <memory>

namespace xtsp
{
  /// @brief Hamiltonian tour
  /// @tparam CostTy (needs to be a signed numeric type)
  template <typename CostTy>
  class Tour
  {
  static_assert(std::is_arithmetic<CostTy>::value && std::is_signed<CostTy>::value);
  public:
    /// @brief Initialize a VALID Hamiltonian tour
    ///         visiting vertices 0, 1, ..., N-1
    ///         based on the permutation \p sequence .
    /// @param sequence a permutation of length N
    /// @param numVertices the variable N  
    Tour(const std::vector<size_t>& sequence, size_t numberVertices);

    // static Tour<CostTy> read_tsplib(std::string_view filepath);
    // virtual void save_tsplib(std::string_view filepath) const;

    /// @brief number of vertex visits in the tour
    virtual size_t size() const final
    {
      return m_seq.size();
    }
  
    /*  not affecting the tour cost */

    /// @brief lower-level access
    /// @param tourRank must be in {0, 1, ..., this->size()-1}
    /// @return the index of the vertex at the requested tour position
    /// @see Tour<CostTy>::getVertex
    virtual size_t getVertex_(size_t tourRank) const final
    {
      return m_seq[tourRank];
    }

    /// @brief higher-level access 
    ///        where we treat the tour as a circular buffer
    ///        (quite useful when implementing sth like 2-opt)
    /// @param tourRank can be any non-negative number,
    /// @return the index of the vertex at the requested tour position
    virtual size_t getVertex(size_t tourRank) const final
    {
      return m_seq[tourRank%size()];
    }

    virtual const CostTy getCost() const final
    {
      return m_cost;
    }

    /// @throw std::invalid_argument
    void assertHamiltonian(size_t numVertices) const;

    /// useful when we need to sort solutions by cost, e.g., in Genetic Algorithms
    bool operator< (const Tour<CostTy>& rhs) const
    {
      return this->m_cost < rhs.m_cost;
    }

    // virtual void reflect();
    // virtual void rightShiftRotate(int amount);
    
    /// ===========================================
    ///  operations that can affect the tour cost
    /// ===========================================

    /// @brief get a mutable reference to the cached cost value.
    ///  Use cases: setting/updating the cost after each 
    ///  tour improvement.  
    virtual CostTy& getCostMutableRef() final
    {
      return m_cost;
    }

    /// @brief get a mutable reference to the sequence.
    /// 
    /// Caller must ensure the tour remains valid.
    /// We provide this function simply for performance reason:
    /// guaranteeing no memory allocation if it's not needed.
    /// 
    /// @see \p xtsp::algo::GtspClusterOptimizer<CostTy>::improve
    virtual std::vector<size_t>& getSeqMutableRef__() final
    {
      return m_seq;
    }

    /// @brief evaluate the cost of a Hamiltonian/generalized tour over the \p graph
    /// @pre Caller shall ensure \p this tour is VALID.
    /// @return immediately retrieve the tour cost 
    ///         (alternatively, you can retrieve the cached value later)
    /// @sa assertHamiltonian
    /// @sa GeneralizedTour<CostTy>::assertGeneralized
    virtual CostTy evalTour(const AbstractCompGraph<CostTy>& graph) final
    {
      // the loop-closing edge
      CostTy cost = graph.getEdgeCost(m_seq[size()-1], m_seq[0]);
      // the rest of the edges
      for (size_t i = 0; i < this->size()-1; ++i)
      {
        cost += graph.getEdgeCost(m_seq[i], m_seq[i+1]);
      }
      m_cost = cost;
      return cost;
    }


    /// @todo update the signature
    /// @todo update the cost
    // virtual void applyTwoOpt(void);

    
  protected:
    // sequence of vertices, m_seq is the "cut" vertex, 
    // which is NOT duplicated at the end of the sequence.
    std::vector<size_t> m_seq;
    // tour cost (we allow it to be "un-initialized")
    CostTy m_cost = std::numeric_limits<CostTy>::max();

    // This constructor only ensures the tour doesn't visit any vertex more than once.
    // Use case: GeneralizedTour's construction
    Tour(const std::vector<size_t>& sequence);

  public: 
    // to remove Wdelete-non-abstract-non-virtual-dtor warning
    virtual ~Tour() = default;
  };

  template <typename CostTy>
  class GeneralizedTour : public Tour<CostTy>
  {
  public:
    /// @brief initialize a VALID generalized tour 
    ///        (over a complete graph implied by \p clustering .)
    /// @throw std::invalid_argument if \p clustering is incompatible 
    ///        with \p vertexSequence .
    GeneralizedTour(
      const std::vector<size_t>& vertexSequence, 
      const std::shared_ptr<Clustering> clustering);

    static GeneralizedTour read_tsplib(
      std::string_view filepath, 
      AbstractCompGraph<CostTy>& g);
    void save_gtsplib(std::string_view filepath) const;
    /// @todo: Karapetyan's format

    size_t numClusters() const;
    size_t numVertices() const;
    
    const std::shared_ptr<Clustering>& getClusteringInfo() const
    {
      return m_clusterInfo;
    }

    void updateCachedClusterSeq();

    /// @brief reverse look-up
    size_t findClusterRankById(size_t clusterId) const;

    size_t getClusterIdByRank(size_t rank) const
    {
      return m_cache_clusterSeq[rank];
    }

    size_t getVertexByClusterId(size_t clusterId) const;

  protected:
    const std::shared_ptr<Clustering> m_clusterInfo;

    /// Lookup-table from rank position to cluster ID.
    /// It shall be updated as long as the object is initialized
    /// and immediately after the cluster sequence is ever changed.
    std::vector<size_t> m_cache_clusterSeq;
  };

} // namespace xtsp

#endif