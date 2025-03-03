#ifndef __XTSP_CORE_TSPLIB_IO_H__
#define __XTSP_CORE_TSPLIB_IO_H__

#include <fstream>
#include <string_view>
#include <Eigen/Core>

#include "xtsp/core/complete_graph.h"

namespace xtsp
{
  enum TsplibFileType
  {
    kTSP,
    kATSP,
    kGTSP,
    kAGTSP,
    kTour,
  };
  // only for problem instances (not for tours)
  enum TsplibEdgeWeightType
  {
    kExplict,
    kEUC_2D, // Euclidean 2D
    kEUC_3D, // Euclidean 3D
    kMAN_2D, // manhatten/L1 norm in 2D
    kMAN_3D, // manhatten/L1 norm in 3D
  };

  /// @throw if unrecognized/unsupported
  enum TsplibFileType tsplibFileTypeFromString(const std::string& val);

  /// @throw if unrecognized/unsupported
  enum TsplibEdgeWeightType tsplibEdgeWeightTypeFromString(const std::string& val);

  /**
   * @brief Parser for files in the TSPLIB format.
   *
   * (it will automatically close the file after the parser is deleted)
   * 
   * The file can be a problem or a tour.
   *
   * In a typical usage, you will use the parser
   * to declare your expectations on the file,
   * in particular, the **order** of data.
   * E.g., For a TSP instance, the file shall
   *
   *   1. "Single-line specifications"
   *      * name
   *      * comment (optional)
   *      * then how many vertices
   *      * EDGE_WEIGHT_TYPE
   *      * ...
   * 
   *   2. data sections
   *      * NODE_COORDINATE_SECTION
   *      * maybe more data section(s)
   *   3. EOF (compulsory?)
   * 
   * Example:
   *
   * The file abc.tsp contains
   * ```
   * NAME : pr342
   * COMMENT: whatever
   * COMMENT: more comment
   * DIMENSION : 342
   * ...
   * NODE_COORDINATE_SECTION
   * 1 0 0
   * 2 3 0
   * ...
   * 342 98 56
   * EOF
   * ```
   * We can use the parser to get the data
   * ```c++
   * TsplibParser parser.("abc.tsp");
   * auto tourName = parser.seekLineAsString("NAME");
   * auto N = parser.seekLineAsInt("DIMENSION");
   * // ... 
   * // suppose we already know it's EUC_2D
   * Eigen::Matrix<int, -1, -1> xy = parser.seekSection(
   *      "NODE_COORDINATE_SECTION", N, true);
   * // ...
   */
  class TsplibParser
  {
  public:
    TsplibParser(const std::string &fpath);

    /**
     * @brief proceed line-by-line until the requested integer is found
     *
     * The requested line should separate the field name and value by a colon.
     * Spaces in the line won't affect the outcome.
     *
     * @tparam T the datatype of the data in question
     * @param fieldName the name of data
     * @return the data in that line
     */
    int seekLineAsInt(const std::string &fieldName);

    /**
     * @brief proceed line-by-line until the requested string is found
     * @see seekLineAsInt
     */
    std::string seekLineAsString(const std::string &fieldName);

    /**
     * @brief seek a floating-point data section
     * 
     * @param numRows the expected number of rows
     * @param numPayloadCols the expected number of payload columns
     * @param hasEnumerationCol used for counting lines in a section 
     *           (not considered part of the payload)
     */
    Eigen::Matrix<float, -1, -1> seekSectionAsFloat(
        const std::string &sectionName,
        size_t numRows, size_t numPayloadCols,
        bool hasEnumerationCol);
    
    /// @param numRows the expected number of clusters 
    /// @param numVertices used to check if all entries are valid
    std::shared_ptr<Clustering> seekGtspSetSection(
        size_t numRows, size_t numVertices);

    void expectReachedEof();

  protected:
    std::ifstream m_fs;
  };
}


#endif