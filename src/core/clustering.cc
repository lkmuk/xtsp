#include "xtsp/core/clustering.h"

#include <algorithm> // std::find
#include <exception>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>
// #include <spdlog/fmt/bundled/ranges.h> // for fmt.join (in debug log)

namespace xtsp
{
  Clustering::Clustering(size_t numVertices, const std::vector<std::vector<size_t>>& membership) 
    : m_n(numVertices), m_c2v(membership) 
  {
    if (numVertices == 0)
    {
      throw std::invalid_argument("It makes no sense to partition an empty graph");
    }
    postInit();
  };

  // Clustering::Clustering(size_t numVertices) 
  //   : m_n(numVertices)
  // {
  //   assert(numVertices > 0);
  //   std::vector <size_t> dummy;
  //   dummy.reserve(numVertices);
  //   for (size_t i = 0; i < numVertices; ++i)
  //     dummy.emplace_back(i);
  //   m_c2v.emplace_back(dummy);
  //   postInit();
  // }

  void Clustering::postInit()
  {
    // we need random access so we need to initialize it now
    size_t uninitializedVal = std::numeric_limits<size_t>::max();
    m_v2c = std::vector<size_t> (this->m_n, uninitializedVal);

    // ii = cluster index
    for (size_t ii = 0; ii < this->numClusters(); ++ii)
    {
      const auto& thisCluster = m_c2v[ii];
      if (thisCluster.empty())
      {
        std::string errMsg = spdlog::fmt_lib::format(
            "Invalid clustering because cluster {:d} is empty.", ii);
        throw std::invalid_argument(errMsg);
      }

      for (const size_t globalVertexId : thisCluster)
      {
        if (globalVertexId >= m_n)
        {
          std::string errMsg = spdlog::fmt_lib::format(
            "Invalid clustering because the maximum vertex index of N-1, "
            "i.e., {:d} is exceeded", 
            globalVertexId);
          throw std::invalid_argument(errMsg);
        }
        size_t& v2cEntry = m_v2c[globalVertexId];
        if (v2cEntry != uninitializedVal)
        {
          std::string errMsg = spdlog::fmt_lib::format(
            "Invalid clustering because each vertex shall be assigned to exactly one cluster "
            "but you assign vertex {:d} to both clusters {:d} and {:d}.", 
            globalVertexId, v2cEntry, ii);
          throw std::invalid_argument(errMsg);
        }
        v2cEntry = ii;
        // spdlog::debug("writing {:d}", ii);
        // spdlog::debug("v2c becomes {}", spdlog::fmt_lib::join(m_v2c, " "));
      }
    }
    
    // ensure no vertex is un-assigned to a cluster
    auto iteVertex = std::find(m_v2c.cbegin(), m_v2c.cend(), uninitializedVal);
    if (iteVertex != m_v2c.cend())
    {
      std::string errMsg = spdlog::fmt_lib::format(
        "Invalid clustering because vertex {:d} is not assigned to any cluster.", 
        std::distance(m_v2c.cbegin(), iteVertex));
      throw std::invalid_argument(errMsg);
    }
  }

  size_t Clustering::getClusterSize(size_t clusterId) const
  {
    return m_c2v[clusterId].size();
  }

  size_t Clustering::evalWhichHasTheLeastVertices() const
  {
    size_t minSoFar = std::numeric_limits<size_t>::max();
    for (size_t i = 0; i < numClusters(); ++i)
    {
      size_t sizeOfThisCluster = getClusterSize(i);
      // early exit
      if (sizeOfThisCluster == 1)
        return i;
      
      if (sizeOfThisCluster < minSoFar)
        minSoFar = sizeOfThisCluster;
    }
    return minSoFar;
  }

  Clustering Clustering::cumsum(const std::vector<size_t>& clustersSizes)
  {
    if (clustersSizes.empty())
      throw std::invalid_argument("clustersSizes should be non-empty");
    
    std::vector<std::vector<size_t>> membership;
    // `head` is the total number of vertices processed so far.
    // it is increment in bulk for each cluster.
    size_t head = 0;
    for (const size_t sizeThisCluster : clustersSizes)
    {
      if (sizeThisCluster == 0)
        throw std::invalid_argument(
          "a cluster is declared as empty, which is not allowed here");
      // member(s) of this cluster
      std::vector<size_t> members;
      for (size_t i = 0; i < sizeThisCluster; ++i)
      {
        members.emplace_back(i+head);
      }
      head = members.back()+1; // update for the next generation
      membership.emplace_back(members);
    }
    
    return Clustering(head, std::move(membership));
  }
} // namespace xtsp