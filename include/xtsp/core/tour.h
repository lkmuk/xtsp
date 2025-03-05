#ifndef __XTSP_CORE_TOUR_H__
#define __XTSP_CORE_TOUR_H__

#include "xtsp/core/complete_graph.h"

#include <string_view>
#include <limits>
#include <type_traits>
#include <memory>

namespace xtsp
{
  /** 
   * Abstract base class for representing a tour
   * 
   * Let seq[i] represents the vertex ID at rank i,
   * i.e., at the (i+1)-th tour position.
   * 
   * Picturally (note that it's directed)
   * ```
   *    +--> seq[0] -> ... -> seq[len(seq)-1] --+
   *    |                                       |
   *    +----------------<----------------------+
   * ```
   * 
   * Some assumptions / non-assumption
   * 
   *  1. It is tour by construction, and 
   *     each modifier must ensure this remains true
   *     when the modication completes.
   * 
   *  2. each vertex will be numbered contiguously
   *     in 0-based indexing:
   * 
   *     from 0,..., \p maxSize() - 1
   * 
   *  2. We allow a tour to be "incomplete" or "partial", 
   *     i.e., \p size() <= \p maxSize() .
   *     In any case, you must specify \p maxSize()
   *     upfront, which you typically can do.
   *  
   * 3. The tour can also revisit some vertex multiple times.
   * 
   */
  struct AbstractTour
  {
    // number of visits in the tour
    virtual size_t size() const = 0;
    // maximum number of visits in the tour
    virtual size_t maxSize() const = 0;

    // the Id of the vertex after vertex v along the tour
    ///
    /// call this only if v is already in the tour
    virtual size_t next(size_t v) const = 0;

    /// the Id of the vertex before vertex v along the tour 
    ///  (Note:  most conventions prefer sticking to "next")
    /// Nonetheless, having both \p prev(v) and \p next(v)
    /// can be quite useful for reinsertion improvement
    ///
    /// call this only if v is already in the tour
    // virtual size_t prev(size_t v) const = 0;

    virtual bool isHamiltonian() const;
    
    virtual bool hasNoRevisit() const;
    virtual bool allElementsAreValid() const;

    // following the implied tour direction
    virtual bool isOneStepAhead(size_t vStart, size_t vGoal) const;
    virtual bool isTwoStepsAhead(size_t vStart, size_t vGoal) const;
    virtual bool isTwoPlusStepsAhead(size_t vStart, size_t vGoal) const;
    virtual size_t evalNumStepsAhead(size_t vStart, size_t vGoal) const;
    
    virtual size_t getDepotId() const = 0;

    /// @brief Save the tour into the TSPLIB format.
    ///
    /// In the process, every vertex ID will be increment by 1 
    /// due to TSPLIB's 1-based indexing scheme.
    ///
    /// What about loading a tour file? 
    /// It's not handled in this base class.
    /// Instead, use PermTour::readTsplib (the "canonical representation).
    /// If needed, convert the resultant PermTour object to another representation. 
    ///
    /// @param fpath output file path (the file extension is upto you)
    /// @param name the NAME of the tour (written in the file)
    /// @todo support COMMENT
    virtual void saveTsplib(const std::string &fpath, const std::string &name) const final;

    ///////////////////////////////
    //   modifiers
    ///////////////////////////////

    /// @brief make vertex \p v rank 0
    /// @todo 
    // virtual void setDepot(size_t v) = 0;

    /**
     * Replace the edges AB and CD in the tour with AC and BD
     * 
     * ```
     *     +... -> A -> B -> ... -> C -> D -> ...+
     *     :                                     :
     *     +.....................................+
     * ```
     * If it reduces the tour cost, it's called a 2-opt move.
     * This operation doesn't care whether it's 2-opt and
     * finding a 2-opt move is outside the scope of this class.
     * The only requirements here are 
     * 
     *   * C must be at least two steps away from A, and
     *   * A must be at least two steps away from C.
     * 
     *   In other words, B != C and A != D.
     * 
     * That said, this operation is mostly meant 
     * for symmetric TSPs and symmetric Generalized TSPs.
     * 
     * Note that this operation requires flipping either 
     * segment BC or AD. If you really want to stipulate 
     * flipping BC (more relevant for ATSP), then 
     * set \p strict as true. Otherwise it's up to the 
     * implementation which one to flip.
     * 
     * @param vA ID of vertex A
     * @param vC ID of vertex C 
     * @param strict See the note above.
     */
    virtual void exchangeTwoEdges(size_t vA, size_t vC, bool strict = false) = 0;

    // virtual void rotateRShift(int v) = 0;
    // virtual void reflectAboutDepot() = 0;
  protected:
    // see Wdelete-non-abstract-non-virtual-dtor
    virtual ~AbstractTour() = default;
  };


  template <typename CostTy>
  CostTy evalTour(const AbstractTour& tour, const AbstractCompGraph<CostTy>& g)
  {
    size_t vHead = tour.getDepotId();
    CostTy sum = 0;
    for (size_t rank = 0; rank < tour.size(); ++rank)
    {
      sum += g.getEdgeCost(vHead, tour.next(vHead));
      vHead = tour.next(vHead);
    }
    return sum;
  }

  /// @brief Permutation representation of a no-revisit tour
  /// 
  ///  where data[i] = vertex ID of at tour position/rank i and
  ///
  ///  Rank index wraps around [0, tour_size - 1].
  ///  
  ///  This rank-based representation is well-suited for book-keeping,
  ///  crossover operators in Genetic Algorithms 
  ///  
  class PermTour : public AbstractTour
  {
  public:
    /// @brief Initialize a tour that has no revisits
    ///         based on the permutation \p sequence .
    /// @param sequence a permutation of length N
    /// @param maxNumVertices the variable N 
    /// @param checks check no revisit and no invalid elements in \p sequence .
    ///               Note if all checks are pass 
    ///               AND sequence.size() == maxNumVertices, then the tour is Hamiltonian
    PermTour(const std::vector<size_t> &sequence, int maxNumVertices = -1, bool checks = true);

    virtual size_t size() const override final
    {
      return m_seq.size();
    }

    virtual size_t maxSize() const override final
    {
      return m_N;
    }

    // Current implementation O(N) access
    virtual size_t next(size_t v) const override;

    virtual size_t getDepotId() const override;

    // O(1) access
    // ID of the vertex after the requested vertex (in rank representation)
    size_t nextByRank(size_t rank) const;

    /// @brief O(n) reverse look-up where n = this->size()
    /// @param vertexId
    /// @param return -1 if not found, otherwise the tour rank
    size_t getRank_(size_t vertexId) const;

    /// @brief lower-level access
    /// @param tourRank must be in {0, 1, ..., this->size()-1}
    /// @return the index of the vertex at the requested tour position
    /// @see Tour<CostTy>::getVertex
    size_t getVertex_(size_t tourRank) const
    {
      return m_seq[tourRank];
    }

    /// @brief higher-level access 
    ///        where we treat the tour as a circular buffer
    ///        (quite useful when implementing sth like 2-opt)
    /// @param tourRank can be any non-negative number,
    /// @return the index of the vertex at the requested tour position
    size_t getVertex(size_t tourRank) const
    {
      return m_seq[tourRank%size()];
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

    /**
     * Swap edges AB and CD in a tour with AC and BD
     *  
     * @param rankA the rank of vertex A  We require 0 <= \p rankA < N
     * @param rankC the rank of vertex C. We require \p rankC >= \p rankA + 2. \p rankC can be arbitraily large (we will wrap it).
     */
    virtual void exchangeTwoEdges_rankBased(size_t rankA, size_t rankC, bool strict = false);
    /**
     * If you are committed to this permuatation representation, 
     * you might want to use exchangeTwoEdges_rankBased instead 
     * @see exchangeTwoEdges_rankBased
     */
    virtual void exchangeTwoEdges(size_t vA, size_t vC, bool strict = false) override;


    
  protected:
    /// sequence of vertices, \p m_seq[0] is the "cut" vertex, 
    /// which is NOT duplicated at the end of the sequence.
    std::vector<size_t> m_seq;
    // max. number of entries in the tour
    const size_t m_N; 

    friend class GeneralizedTour;
  };

  /// Performance-oriented overload 
  /// 
  /// For array representation, 
  /// the general implementation for abstract tour is a bit inefficient.
  /// We use this for better performance (avoid reverse lookup)
  /// even skipping the modulo operation in \p PermTour::getVertex()
  template <typename CostTy>
  CostTy evalTour(const PermTour &tour, const AbstractCompGraph<CostTy>& g)
  {
    CostTy sum = g.getEdgeCost(tour.getVertex_(tour.size()-1), tour.getVertex_(0));
    for (size_t rank = 1; rank < tour.size(); ++rank)
      sum += g.getEdgeCost(tour.getVertex_(rank-1), tour.getVertex_(rank));
    return sum;
  }

  /// @todo use the unified API (e.g., make it a template?)
  /// @todo What about thread-safety???
  class GeneralizedTour
  {
  public:
    /// @brief initialize a VALID generalized tour 
    ///        (over a complete graph implied by \p clustering .)
    /// @throw std::invalid_argument if \p clustering is incompatible 
    ///        with \p tour .
    ///
    GeneralizedTour(
      std::shared_ptr<PermTour> tour, 
      const std::shared_ptr<Clustering> clustering,
      bool check = true);

    // sometimes this API is more friendly than the native constructor
    static GeneralizedTour fromPermutation(
      const std::vector<size_t>& tour, 
      const std::shared_ptr<Clustering> clustering,
      bool check = true);

    size_t numClusters() const;
    size_t numVertices() const;
    
    const std::shared_ptr<Clustering>& getClusteringInfo() const
    {
      return m_clusterInfo;
    }

    const std::shared_ptr<PermTour> getTour() const;
    const std::shared_ptr<PermTour> getSuperTour() const;
    std::shared_ptr<PermTour> getTourMutable__();
    std::shared_ptr<PermTour> getSuperTourMutable__();

    void updateCachedSuperTour();

    /// @todo get rid of this API
    size_t getClusterIdByRank(size_t rank) const
    {
      return m_cache_supTour->getVertex(rank);
    }

    /// @brief reverse look-up
    /// @todo needed?
    size_t findClusterRankById(size_t clusterId) const;

    size_t getVertexByClusterId(size_t clusterId) const;

  protected:
    std::shared_ptr<PermTour> m_tour;
    const std::shared_ptr<Clustering> m_clusterInfo;

    /// Lookup-table from rank position to cluster ID.
    /// It shall be updated as long as the object is initialized
    /// and immediately after the cluster sequence is ever changed.
    std::shared_ptr<PermTour> m_cache_supTour = nullptr; // initialized in updateCachedSuperTour

  };

} // namespace xtsp

#endif