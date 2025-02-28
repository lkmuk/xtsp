#ifndef __XTSP_CORE_CLUSTERING_H__
#define __XTSP_CORE_CLUSTERING_H__

#include <vector>
#include <stddef.h> // size_t

namespace xtsp
{

  ///@brief useful for Generalized TSP, 
  ///@note In the simple case, just take a single cluster for all graph vertices
  class Clustering
  {
  public:
    /// @brief create a non-trivial clustering
    /// @param numVertices 
    /// @param membership membership[i] is the set of vertices associated to cluster i
    Clustering(size_t numVertices, const std::vector<std::vector<size_t>>& membership);

    /// @brief create a trivial clustering (a single cluster for all vertices)
    /// @param numVertices 
    // Clustering(size_t numVertices);

    const size_t numVertices() const
    {
      return m_v2c.size();
    }

    const size_t numClusters() const
    {
      return m_c2v.size();
    }

    /// @brief return the set of indices of vertices belonging to a cluster
    const std::vector<size_t>& getMembers(size_t clusterId) const
    {
      return m_c2v[clusterId];
    }

    /// @brief reverse lookup from global vertex ID to cluster ID
    size_t getClusterId(size_t vertex) const
    {
      return m_v2c[vertex];
    }


  protected:
    // total number of vertices (among all clusters)
    size_t m_n;
    // m_c2v[i] = the list of all its vertices' indices of cluster i
    std::vector<std::vector<size_t>> m_c2v;
    // pre-computed reverse lookup from vertex ID to cluster ID
    std::vector<size_t> m_v2c;

    /// @brief populate the reverse LUT `m_v2c` while checking if the "clustering" partitions 0, ..., N-1 
    /// @note meant to be called during the construction
    void postInit();
  };
} // namespace

#endif