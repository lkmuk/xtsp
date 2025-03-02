#include "xtsp/core/complete_graph.h"

#include <Eigen/Dense> // array ops
#include <exception>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>

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
  const Eigen::Matrix<CostTy, -1, -1>& ImplicitCompleteGraph<CostTy>::getXy() const
  {
    return m_xy;
  }

  template <typename CostTy>
  size_t ImplicitCompleteGraph<CostTy>::nDim() const
  {
    return m_xy.cols();
  }

  // template <typename CostTy>
  // AbstractCompGraph<CostTy> loadTsplibProblem(
  //   std::string_view fpath, bool isGeneralized)
  // {

  // }

  // explicit instantiation
  template class CompleteGraph<float>;
  template class CompleteGraph<int>;
  template class ImplicitCompleteGraph<float>;
} // namespace xtsp