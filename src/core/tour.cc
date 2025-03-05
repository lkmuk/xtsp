#include "xtsp/core/tour.h"
#include "xtsp/core/tour_alternatives.h"

#include "xtsp/core/utils.h"
#include "../toolbox/ring_ops.h"
#include "../toolbox/work_buffer.h"

#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>

namespace xtsp
{
  void AbstractTour::print (std::ostream &os) const
  {
    auto vHead = getDepotId();
    for (size_t rank = 0; rank < size(); ++ rank)
    {
      os << vHead << "-";
      vHead = next(vHead);
    }
  }


  bool AbstractTour::isOneStepAhead(size_t vStart, size_t vGoal) const
  {
    return next(vStart) == vGoal;
  }

  bool AbstractTour::isTwoStepsAhead(size_t vStart, size_t vGoal) const
  {
    return next(next(vStart)) == vGoal;
  }

  bool AbstractTour::isTwoPlusStepsAhead(size_t vStart, size_t vGoal) const
  {
    return (!isOneStepAhead(vStart, vGoal)) && (vStart != vGoal);
  }

  size_t AbstractTour::evalNumStepsAhead(size_t vStart, size_t vGoal) const
  {
    size_t count = 0;
    size_t head = vStart;
    while (head != vGoal)
    {
      ++count;
      head = next(head);
      assert(count < size()); // to catch bug
    }
    return count;
  }

  bool AbstractTour::isHamiltonian() const
  {
    // Note: a partial tour is not Hamiltonian
    return (size() < maxSize()) ? false : (hasNoRevisit() && allElementsAreValid());
  }

  bool AbstractTour::allElementsAreValid() const
  {
    size_t vHead = getDepotId();
    for (size_t rank = 0; rank < size(); ++rank)
    {
      if (vHead >= maxSize())
        return false;
      vHead = next(vHead);
    }
    return true;
  }


  bool AbstractTour::hasNoRevisit() const
  {
    // While we can use std::set, the
    // O(1) random access offered by WorkBuffer is 
    // much better performance wise.
    //
    // all entries are initially marked as false,
    // i.e., not found.
    internal::WorkBuffer buf(maxSize(), false);

    size_t vHead = getDepotId();
    for (size_t rank = 0; rank < size(); ++rank)
    {
      // recall TODO == not found in tour
      //      DONE == found in tour
      if (buf.isDone(vHead))
        return false; 
      buf.markDone(vHead);
      vHead = next(vHead);
    }
    return true;
  }

  ////////////////////////////////////////
  ///   PermTour : ctor & key interfaces
  ////////////////////////////////////////

  PermTour::PermTour(const std::vector<size_t> &sequence, int maxNumVertices, bool checks)
    : m_seq(sequence)
    , m_N((maxNumVertices < 0) ? sequence.size() : maxNumVertices)
  {
    if (sequence.size() > m_N)
      throw std::invalid_argument(
        "[constructing a PermTour] maxNumVertices should >= sequence.size()");
    if (checks)
    {
      xtsp::utils::assertNoDuplicate(sequence, "tour", "city");
      /// @todo unit test coverage
      xtsp::utils::assertAllValid(m_N, sequence, "tour", "city ID");
    }
  }

  size_t PermTour::next(size_t v) const
  {
    return nextByRank(getRank_(v));
  }

  size_t PermTour::getDepotId() const
  {
    return m_seq[0];
  }

  void PermTour::exchangeTwoEdges_rankBased(size_t rankA, size_t rankC, bool strict)
  {
    if (rankA + 1 >= rankC)
      throw std::invalid_argument("expect rankA + 1 < rankC");
    size_t rankB = rankA + 1;
    if (strict)
      xtsp::internal::reverseRingSegment_strict(m_seq, rankB, rankC);
    else
      xtsp::internal::reverseRingSegment_smart(m_seq, rankB, rankC);
  }
  void PermTour::exchangeTwoEdges(size_t vA, size_t vC, bool strict)
  {
    if (vA == vC)
    {
      SPDLOG_WARN("Ignoring a no-op two-edge-exchange request");
    }
    // some overhead as you can see
    size_t rankA = this->getRank_(vA);
    size_t rankC = this->getRank_(vC);
    // some quirks
    if (rankC < rankA)
      rankC += size();
    // SPDLOG_ERROR("rankA = {:d}, rankC = {:d}", rankA, rankC);
    exchangeTwoEdges_rankBased(rankA, rankC, strict);
  }

  ////////////////////////////////////////////////////
  ///   Some native getters of PermTour
  ////////////////////////////////////////////////////

  size_t PermTour::getRank_(size_t vertexId) const
  {
    auto ite = std::find(m_seq.cbegin(), m_seq.cend(), vertexId);
    if (ite == m_seq.cend())
      throw ("Cannot find the requested vertex ID in the tour");
    return std::distance(m_seq.cbegin(), ite);
  }


  size_t PermTour::nextByRank(size_t rank) const
  {
    return m_seq[(rank+1)%size()];
  }

  ////////////////////////////////////////////
  ///  Adjacency table representation
  ////////////////////////////////////////////
  AdjTabTour::AdjTabTour(const std::vector<size_t>& perm, int maxNumVertices, bool checks)
    : m_N((maxNumVertices < 0) ? perm.size() : maxNumVertices)
  {
    if (perm.size() == 0)
      throw std::invalid_argument(
        "[constructing a AdjTabTour] the input permutation vector is empty");
    m_head = perm[0];
    if (perm.size() > m_N)
      throw std::invalid_argument(
        "[constructing a PermTour] maxNumVertices should >= sequence.size()");
    if (checks)
    {
      xtsp::utils::assertNoDuplicate(perm, "tour", "city");
      xtsp::utils::assertAllValid(m_N, perm, "tour", "city ID");
    }

    m_cache_tourSize = perm.size();
    m_dat.reserve(m_N);
    // a simple implementation that ensures that 
    // those not yet in the tour also get properly initialized
    for (size_t v = 0; v < m_N; ++v)
    {
      m_dat.emplace_back(v);
    }
    for (size_t rank = 0; rank < m_cache_tourSize; ++rank)
    {
      size_t vertexId = perm[rank];
      auto rankPrevWrapped = (rank + m_cache_tourSize - 1)%m_cache_tourSize;
      auto rankNextWrapped = (rank + m_cache_tourSize + 1)%m_cache_tourSize;
      m_dat[vertexId].prev = m_dat.data() + perm[rankPrevWrapped];
      m_dat[vertexId].next = m_dat.data() + perm[rankNextWrapped];
    }
    
  }
  size_t AdjTabTour::size() const
  {
    return m_cache_tourSize;
  }
  size_t AdjTabTour::maxSize() const
  {
    return m_cache_tourSize;
  }

  size_t AdjTabTour::next(size_t vertex) const
  {
    return m_dat[vertex].next->id;
  }

  size_t AdjTabTour::getDepotId() const
  {
    return m_head;
  }

  void AdjTabTour::exchangeTwoEdges(
      size_t vA, size_t vC, bool strict)
  {
    // (only in Debug build)
    /// ensure  everything is valid 
    assert(isTwoPlusStepsAhead(vA, vC));
    assert(isTwoPlusStepsAhead(vC, vA));

    const size_t vB = next(vA);
    const size_t vD = next(vC);

    /// currently we always reverse CD
    /// i.e., ignoring the strict-flag

    // if the segment CD has some node(s) in between, then ...
    // modify the underlying data while traversing along 
    // the (original) tour direction.
    size_t vHead = next(vC);
    while (vHead != vD) 
    {
      std::swap(m_dat[vHead].prev, m_dat[vHead].next);
      vHead = m_dat[vHead].prev->id; // i.e., next before the flip
    }

    m_dat[vA].next = &m_dat[vC]; // was vB
    m_dat[vD].prev = &m_dat[vB]; // was vC

    // change vB, ..., vC
    m_dat[vB].prev = m_dat[vB].next;
    m_dat[vB].next = &m_dat[vD];

    m_dat[vC].next = m_dat[vC].prev;
    m_dat[vC].prev = &m_dat[vA];
  }

  


  ////////////////////////////////////////////
  ///  GeneralizedTour
  ////////////////////////////////////////////

  GeneralizedTour GeneralizedTour::fromPermutation(
    const std::vector<size_t>& tour, 
    const std::shared_ptr<Clustering> clustering,
    bool check)
  {
    return GeneralizedTour(
      std::make_shared<PermTour>(tour, clustering->numVertices()), 
      clustering, 
      check);
  }

  GeneralizedTour::GeneralizedTour(
    std::shared_ptr<PermTour> tour, 
    const std::shared_ptr<Clustering> clustering,
    bool check) 
    : m_tour(tour), m_clusterInfo(clustering)
  {
    updateCachedSuperTour();
    if (check)
    {
      xtsp::utils::assertNoDuplicate(m_cache_supTour->m_seq, "generalized tour", "cluster");
    }
  }

  void GeneralizedTour::updateCachedSuperTour()
  {
    if (m_cache_supTour == nullptr)
    {
      // dummy allocation
      m_cache_supTour = std::make_shared<PermTour>(
        std::vector<size_t>({0}), m_clusterInfo->numClusters(), false);
      
      m_cache_supTour->m_seq.reserve(m_clusterInfo->numClusters());
    }

    m_cache_supTour->m_seq.clear();

    // SPDLOG_ERROR("m_tour ptr {}", (void *)m_tour.get());

    /// @todo use the unified API
    for (const size_t vertexId : this->m_tour->getSeqMutableRef__())
    {
      size_t whichCluster = m_clusterInfo->getClusterId(vertexId);
      m_cache_supTour->m_seq.emplace_back(whichCluster);
    }
  }

  size_t GeneralizedTour::numClusters() const
  {
    return m_clusterInfo->numClusters();
  }

  size_t GeneralizedTour::numVertices() const
  {
    return m_clusterInfo->numVertices();
  }

  const std::shared_ptr<PermTour> GeneralizedTour::getTour() const
  {
    return m_tour;
  }
  const std::shared_ptr<PermTour> GeneralizedTour::getSuperTour() const
  {
    return m_cache_supTour;
  }
  std::shared_ptr<PermTour> GeneralizedTour::getTourMutable__()
  {
    return m_tour;
  }
  std::shared_ptr<PermTour> GeneralizedTour::getSuperTourMutable__()
  {
    return m_cache_supTour;
  }

  size_t GeneralizedTour::findClusterRankById(size_t clusterId) const
  {
    if (clusterId >= this->numClusters())
      throw std::out_of_range("the requested cluster ID is too large");
    auto ite = std::find(m_cache_supTour->m_seq.cbegin(), m_cache_supTour->m_seq.cend(), clusterId);
    size_t result = std::distance(m_cache_supTour->m_seq.cbegin(), ite);
    assert(result < this->numClusters());
    return result;
  }

  size_t GeneralizedTour::getVertexByClusterId(size_t clusterId) const
  {
    size_t clusterRank = findClusterRankById(clusterId);
    return this->m_tour->m_seq[clusterRank];
  }

} // namespace xtsp