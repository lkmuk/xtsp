#include "xtsp/core/complete_graph.h"

#include <Eigen/Dense> // array ops
#include <exception>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>

#include "xtsp/core/tsplib_io.h"

namespace xtsp
{
  template <typename CostTy>
  AbstractCompGraph<CostTy>::AbstractCompGraph(const std::shared_ptr<Clustering> clustering)
    : m_clustering(clustering)
  {
  }

  template <typename CostTy>
  size_t AbstractCompGraph<CostTy>::numClusters() const
  {
    return (m_clustering == nullptr) ? 1 : m_clustering->numClusters();
  }

  template <typename CostTy>
  bool AbstractCompGraph<CostTy>::isClustered() const
  {
    return !(m_clustering == nullptr);
  }

  template <typename CostTy>
  CompleteGraph<CostTy>::CompleteGraph(
    bool isSymmetric, 
    const Eigen::Matrix<CostTy, -1, -1>& edgeCosts, 
    const std::shared_ptr<Clustering> clustering)
    : AbstractCompGraph<CostTy>(clustering), 
      m_symmetric(isSymmetric),
      m_mat(edgeCosts)
  {
    if (m_mat.rows() != m_mat.cols())
    {
      throw std::invalid_argument("The edge cost matrix must be square!");
    }

    // https://stackoverflow.com/questions/13403982/eigen-efficient-type-for-dense-symmetric-matrix
    // somehow fail to compile
    // m_mat = (m_symmetric) ? edgeCosts.selfadjointView<Eigen::Lower>() : edgeCosts;
    if (isSymmetric)
    {
      // enforce symmetry 
      // (by copying from the lower triangle to the upper triangle)
      for (size_t i = 0; i < static_cast<size_t>(m_mat.rows()) - 1; ++i)
        for (size_t j = i+1; j < static_cast<size_t>(m_mat.cols()); ++j)
          m_mat(i,j) = m_mat(j,i);
    }

    /// check also if >= 0 ?
    for (size_t i = 0; i < static_cast<size_t>(m_mat.rows()); ++i)
    {
      for (size_t j = 0; j < static_cast<size_t>(m_mat.cols()); ++j)
      {
        if (m_mat(i,j) < static_cast<CostTy>(0))
        {
          spdlog::warn(
            "The cost for edge {:d}->{:d} is negative. "
            "Some algorithms/tools will fail",
            i, j
          );
        }
      }
    }
  }

  template <typename CostTy>
  bool CompleteGraph<CostTy>::isSymmetric() const
  {
      return m_symmetric;
  }

  template <typename CostTy>
  ImplicitCompleteGraph<CostTy>::ImplicitCompleteGraph(
    const Eigen::Matrix<CostTy, -1, -1>& xy, 
    const std::shared_ptr<Clustering> clustering,
    int normTy)
    : AbstractCompGraph<CostTy>(clustering), 
      m_xy(xy),
      m_normType(normTy)
  {
    if (xy.cols() == 0)
      throw std::invalid_argument("xy data cannot be empty (e.g., having no column)");
    if (normTy > 3 || normTy < 0)
      throw std::invalid_argument("norm type must be 0 or 1 or 2.");
    /// @todo perhaps preevaluate the cost for small-instance problems?
    /// @todo precompute K-d tree
  }

  template <typename CostTy>
  ImplicitCompleteGraph<CostTy> ImplicitCompleteGraph<CostTy>::loadFromTsplibFile(
      const std::string& fpath)
  {
    TsplibParser parser(fpath);
    std::string errMsg;
    std::string problemName = parser.seekLineAsString("NAME");
    SPDLOG_INFO("Parsing tsplib file: NAME = {}", problemName);
    std::string tspTypeStr = parser.seekLineAsString("TYPE");
    SPDLOG_INFO("Parsing tsplib file: TYPE = {}", tspTypeStr);
    enum TsplibFileType tspType = tsplibFileTypeFromString(tspTypeStr);
    bool isGeneralized = false;
    switch (tspType)
    {
      case TsplibFileType::kGTSP:
        isGeneralized = true;
        break;
      case TsplibFileType::kTSP:
        isGeneralized = false;
        break;
      default:
        errMsg = fmt::format(
          "TYPE {} is recognized but not compatible here.", tspTypeStr);
        SPDLOG_ERROR(errMsg);
        throw std::invalid_argument(errMsg);
    }

    int numVertices = parser.seekLineAsInt("DIMENSION");
    SPDLOG_INFO("Parsing tsplib file: DIMENSION = {}", numVertices);
    if (numVertices <= 0)
      throw std::invalid_argument("Bad DIMENSION: should be a positive number");

    int numClusters = 0;
    if (isGeneralized)
    {
      numClusters = parser.seekLineAsInt("GTSP_SETS");
      SPDLOG_INFO("Parsing tsplib file: GTSP_SETS = {}", numClusters);
      if (numClusters < 2)
        throw std::invalid_argument("Bad GTSP_SETS: should be at least 2");
    }

    size_t nDim = 0;
    int normType = -100;
    std::string edgeWeightTypeStr = parser.seekLineAsString("EDGE_WEIGHT_TYPE");
    enum TsplibEdgeWeightType edgeWeightType 
      = tsplibEdgeWeightTypeFromString(edgeWeightTypeStr);
    switch (edgeWeightType)
    {
      case TsplibEdgeWeightType::kEUC_2D:
        nDim = 2;
        normType = 2;
        break;
      case TsplibEdgeWeightType::kEUC_3D:
        nDim = 3;
        normType = 2;
        break;
      case TsplibEdgeWeightType::kMAN_2D:
        nDim = 2;
        normType = 1;
        break;
      case TsplibEdgeWeightType::kMAN_3D:
        nDim = 3;
        normType = 1;
        break;
      default:
        errMsg = fmt::format(
          "EDGE_WEIGHT_TYPE {} is recognized but not compatible here.", edgeWeightTypeStr);
        throw std::invalid_argument(errMsg);
    }

    Eigen::Matrix<CostTy,-1,-1> xy = parser.seekSectionAsFloat(
      "NODE_COORD_SECTION", numVertices, nDim, true);
    SPDLOG_INFO("Parsing tsplib file: successfully parsed the NODE_COORD_SECTION.");

    std::shared_ptr<Clustering> clustering = nullptr;
    if (isGeneralized)
      clustering = parser.seekGtspSetSection(numClusters, numVertices);
    
    parser.expectReachedEof();
    return ImplicitCompleteGraph(xy, clustering, normType);
  }

  template <typename CostTy>
  const Eigen::Matrix<CostTy, -1, -1>& ImplicitCompleteGraph<CostTy>::getXy() const
  {
    return m_xy;
  }

  template <typename CostTy>
  size_t ImplicitCompleteGraph<CostTy>::nDim() const
  {
    return m_xy.cols();
  }

  template <typename CostTy>
  CompleteGraph<int> ImplicitCompleteGraph<CostTy>::explicitize(float scale) const
  {
    if (numVertices() > 10000)
      SPDLOG_WARN(
        "You are trying to create a N x N matrix with large N > 10000. "
        "The conversion may fail due to limited memory.");
    if (scale < 1)
      SPDLOG_WARN(
        "you choose scale = {:.f} but you probably want e.g., scale >= 100", scale);

    Eigen::Matrix<int, -1, -1> costMat(numVertices(), numVertices());
    costMat.diagonal().setZero();
    for (int i = 0; i < (int)numVertices(); ++i)
      for (int j = 0; j < i; ++j)
      {
        int e = nint(scale * getEdgeCost(i,j));
        costMat(i,j) = e; // the lower triangle
        costMat(j,i) = e; // the upper triangle
      }
    return CompleteGraph<int>(true, costMat, this->m_clustering);    
  }


  template <typename CostTy>
  ImplicitCompleteGraph<CostTy> ImplicitCompleteGraph<CostTy>::buildClusterMeans() const
  {
    std::string errMsg;
    if (!this->isClustered())
    {
      errMsg = "The graph is not clustered so buildClusterMeans fails";
      SPDLOG_ERROR(errMsg); // to locate the line
      throw std::invalid_argument(errMsg);
    }
    if (m_normType != 2)
    {
      SPDLOG_WARN(
        "Cluster centroids derived from averaging may not be meaningful "
        " when the edge cost isn't L2-norm.");
    }
    Eigen::Matrix<CostTy, -1, -1> mean (this->numClusters(), nDim());
    mean.setZero();
    for (size_t m = 0; m < this->numClusters(); ++m)
    {
      // since a cluster's vertices may be non-contiguous, 
      // we have to first add the vertices, then average the sum
      for (const size_t n : this->m_clustering->getMembers(m))
        mean.row(m) += m_xy.row(n);
      mean.row(m) /= this->m_clustering->getClusterSize(m);
    }

    return ImplicitCompleteGraph<CostTy> (mean, this->m_clustering, this->m_normType);
  }


  // explicit instantiation
  template class CompleteGraph<float>;
  template class CompleteGraph<int>;
  template class ImplicitCompleteGraph<float>;
} // namespace xtsp