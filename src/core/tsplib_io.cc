#include "xtsp/core/tsplib_io.h"
#include "xtsp/core/tour.h"

#include "tsplib_io_seek_impl.h"
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/core.h>
#include <spdlog/fmt/bundled/ranges.h>

#include <fstream>

typedef std::string string;

namespace xtsp
{

  enum TsplibFileType tsplibFileTypeFromString(const std::string& val)
  {
    if (val == "TSP")
      return TsplibFileType::kTSP;
    if (val == "ATSP")
      return TsplibFileType::kATSP;
    if (val == "GTSP")
      return TsplibFileType::kGTSP;
    if (val == "AGTSP")
      return TsplibFileType::kAGTSP;
    if (val == "TOUR")
      return TsplibFileType::kTour;
    std::string errMsg = fmt::format("Unrecognized TSPLIB file type: {}", val);
    throw std::invalid_argument(errMsg);
  }

  enum TsplibEdgeWeightType tsplibEdgeWeightTypeFromString(const std::string& val)
  {
    if (val == "EXPLICIT")
      return TsplibEdgeWeightType::kExplict;
    if (val == "EUC_2D")
      return TsplibEdgeWeightType::kEUC_2D;
    if (val == "EUC_3D")
      return TsplibEdgeWeightType::kEUC_3D;
    if (val == "MAN_2D")
      return TsplibEdgeWeightType::kMAN_2D;
    if (val == "MAN_3D")
      return TsplibEdgeWeightType::kMAN_3D;
    std::string errMsg = fmt::format("Unrecognized TSPLIB EDGE_WEIGHT_TYPE: {}", val);
    throw std::invalid_argument(errMsg);
  }



  TsplibParser::TsplibParser(const string& fpath)
    : m_fs(std::ifstream(fpath))
  {
    if (!m_fs.is_open() || !m_fs.good())
    {
      std::string errMsg = fmt::format("Failed to open the file: {}", fpath);
      SPDLOG_ERROR(errMsg);
      throw std::runtime_error(errMsg);
    }
    SPDLOG_DEBUG("Opening file {}", fpath);
  }

  int TsplibParser::seekLineAsInt(const string& fieldName)
  {
    return std::atoi(internal::seekLine(m_fs, fieldName).c_str());
  }

  string TsplibParser::seekLineAsString(const string& fieldName)
  {
    return internal::seekLine(m_fs, fieldName);
  }

  void TsplibParser::expectReachedEof()
  {
    std::string line;
    if ((!std::getline(m_fs, line)) && line != "EOF")
    {
      SPDLOG_WARN(
        "Expecting file to end but it doesn't. "
        "We ignore starting from this line:\n{}",
        line);
    }
  }

  Eigen::Matrix<float, -1, -1> TsplibParser::seekSectionAsFloat(
      const string &sectionName,
      size_t numRows, size_t numPayloadCols,
      bool hasEnumerationCol)
  {
    return internal::seekSection(m_fs, sectionName, numRows, numPayloadCols, hasEnumerationCol);
  }

  enum tsplibGtspSetSectionParseState
  {
    kNewline,
    kScanningData,
  };

  std::shared_ptr<Clustering> TsplibParser::seekGtspSetSection(
        size_t numRows, size_t numVertices)
  {
    SPDLOG_INFO("Parsing GTSP_SET_SECTION");
    using string = std::string;
    string line;
    string errMsg;
    if (numRows < 2)
      throw std::invalid_argument("numRows must be >=2");
    while (std::getline(m_fs, line))
    {
      if (line.find("GTSP_SET_SECTION") == std::string::npos)
        continue;
      else // found the section, time to scan the rest
      {
        std::vector<std::vector<size_t>> memberships(numRows);
        int currentClusterId = 1; // 1-based
        int buf;
        enum tsplibGtspSetSectionParseState state = kNewline;
        while (m_fs >> buf)
        {
          switch (state)
          {
            case kNewline:
              if (buf != currentClusterId)
                SPDLOG_WARN("The row for cluster {:d} (1-based index) begins with an unexpected value {:d}", currentClusterId, buf);
              state = kScanningData;
              break;
            case kScanningData:
              if (buf == -1)
              {
                if (++currentClusterId > (int)numRows)
                {
                  SPDLOG_INFO("Parsing GTSP_SET_SECTION : Processed all rows");                  
                  return std::make_shared<Clustering>(numVertices, memberships);
                }
                int lastClusterId = currentClusterId - 2; // zero-based indexing
                SPDLOG_DEBUG(
                  "Converted to 0-based indexing: clusterId = {:d}, members = [{}]", 
                  lastClusterId, fmt::join(memberships[lastClusterId], " "));
                state = kNewline;
              }
              else if (1 <= buf && buf <= (int)numVertices)
              {
                memberships[(size_t)(currentClusterId-1)].emplace_back((size_t)buf-1);
              }
              else
              {
                errMsg = fmt::format(
                  "Unexpected value {:d} while scanning for cluster {:d} (1-indexing): permitted values are -1 or 1, ..., {:d}", buf, currentClusterId, numVertices);
                throw std::invalid_argument(errMsg);
              }
              break;
          }
        }
        if (currentClusterId <= (int)numRows)
        {
          throw std::invalid_argument("Incomplete GTSP_SET_SECTION");
        }
      }
    }
    errMsg = "Failed to find GTSP_SET_SECTION.";
    throw std::invalid_argument(errMsg);
  }

  void AbstractTour::saveTsplib(const std::string &fpath, const std::string &name) const
  {
    std::ofstream fs (fpath);
    fs << "NAME : " << name << "\n";
    fs << "TYPE : TOUR\n";
    fs << "DIMENSION : " << size() << "\n";
    fs << "TOUR_SECTION\n";
    size_t vHead = getDepotId();
    for (size_t rank = 0; rank < size(); ++rank)
    {
        // converting to 1-based indexing
        fs << vHead + 1 << "\n";
        vHead = next(vHead);
    }
    fs << "-1\n";
    fs << "EOF\n";
    fs.close();
  }
}