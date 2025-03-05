#include "xtsp/core/tour.h"

#include "xtsp/core/utils.h"
#include "../toolbox/ring_ops.h"
#include "../toolbox/work_buffer.h"

#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>

namespace xtsp
{
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
    if (rankA + 1 > rankC)
      throw std::invalid_argument("expect rankA + 1 > rankC");
    if (strict)
      xtsp::internal::reverseRingSegment_strict(m_seq, rankA, rankC);
    else
      xtsp::internal::reverseRingSegment_smart(m_seq, rankA, rankC);
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