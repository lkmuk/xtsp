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

    size_t getClusterSize(size_t clusterId) const;

    /// @brief which cluster has the least vertices?
    /// @return cluster index (on tie-break, return the lowest cluster Id)
    size_t evalWhichHasTheLeastVertices() const;
    
    /**
     * @brief Derive the global vertices and partition from "a distribution"
     * 
     * Merge/flatten local vertex indices (which are implicitly) 
     * into a global index space, so it's like **cumsum** 
     * in NumPy and MATLAB
     * 
     * Each cluster occupies a contiguous global index space, 
     * A cluster with a lower cluster index sits lower in the 
     * global space as well. 
     * 
     * Example: Suppose \p clustersSizes is {3, 1, 2}
     *  The flattened indices will be a rather boring 0,1,2,3,4,5
     *  with the cluster membership {{0,1,2},{3}, {4, 5}}.
     *
     * This is useful when we construct a partitioned graph bottom-up.
     * @param clustersSizes It is like a distribution, 
     *    the m-th element is the number of vertices in cluster m.
     *    It should not be empty
     * @throw std::invalid_argument when there is no cluster, or a cluster is empty.
     */
    static Clustering cumsum(const std::vector<size_t>& clustersSizes); 

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