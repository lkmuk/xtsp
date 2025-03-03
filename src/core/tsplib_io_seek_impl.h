#pragma once

#include <fstream>
#include <Eigen/Core>

namespace xtsp::internal
{
  std::string seekLine(std::ifstream& fs, const std::string& fieldName);

  Eigen::Matrix<float,-1,-1> seekSection(
    std::ifstream& fs, const std::string& sectionName, 
    size_t numRows, size_t numPayloadCols, bool hasEnumerationCol);
}