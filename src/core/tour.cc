#include "xtsp/core/tour.h"

#include "xtsp/core/utils.h"

#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>

namespace xtsp
{

  template <typename CostTy>
  Tour<CostTy>::Tour(const std::vector<size_t>& sequence, size_t numVertices)
    : m_seq(sequence)
  {
    assertHamiltonian(numVertices);
  }

  template <typename CostTy>
  Tour<CostTy>::Tour(const std::vector<size_t>& sequence)
    : m_seq(sequence)
  {
    xtsp::utils::assertNoDuplicate(sequence, "tour", "city");
  }



  template <typename CostTy>
  void Tour<CostTy>::assertHamiltonian(size_t numVertices) const
  {
    /// the implementation is quite similar to Clustering::postInit

    if (numVertices == 0)
    {
      throw std::invalid_argument("numVertices should be >= 0.");
    }
    xtsp::utils::assertIsPermutation(numVertices, m_seq, "tour", "city");
    
  }

  template <typename CostTy>
  GeneralizedTour<CostTy>::GeneralizedTour(
    const std::vector<size_t>& sequence, 
    const std::shared_ptr<Clustering> clustering)
    : Tour<CostTy>(sequence), m_clusterInfo(clustering)
  {
    m_cluster2Vertices = std::vector<size_t>(m_clusterInfo->numClusters(), 0);

    for (size_t rank = 0; rank < this->size(); ++rank)
    {
      size_t vertexId = this->m_seq[rank];
      size_t whichCluster = m_clusterInfo->getClusterId(vertexId);
      m_clusterSeq.emplace_back(whichCluster);
      m_cluster2Vertices[whichCluster] = vertexId;
    }

    xtsp::utils::assertNoDuplicate(m_clusterSeq, "generalized tour", "cluster");
  }

  template <typename CostTy>
  size_t GeneralizedTour<CostTy>::numClusters() const
  {
    return m_clusterInfo->numClusters();
  }

  template <typename CostTy>
  size_t GeneralizedTour<CostTy>::numVertices() const
  {
    return m_clusterInfo->numVertices();
  }

  template <typename CostTy>
  size_t GeneralizedTour<CostTy>::findClusterRankById(size_t clusterId) const
  {
    if (clusterId >= this->numClusters())
      throw std::invalid_argument("the requested cluster ID is too large");
    auto ite = std::find(m_clusterSeq.cbegin(), m_clusterSeq.cend(), clusterId);
    size_t result = std::distance(m_clusterSeq.cbegin(), ite);
    assert(result < this->numClusters());
    return result;
  }

  // explicit template instantiation
  template class Tour<float>;
  template class Tour<int>;
  template class GeneralizedTour<float>;
  template class GeneralizedTour<int>;
} // namespace xtsp