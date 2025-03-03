#include "tsplib_io_seek_impl.h"

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>
#include <exception>

namespace xtsp::internal
{
  // https://stackoverflow.com/questions/1798112/removing-leading-and-trailing-spaces-from-a-string
  inline std::string &ltrim(std::string &s, const char *t = " \t\n\r\f\v")
  {
    s.erase(0, s.find_first_not_of(t));
    return s;
  }
  // trim from right
  inline std::string &rtrim(std::string &s, const char *t = " \t\n\r\f\v")
  {
    s.erase(s.find_last_not_of(t) + 1);
    return s;
  }
  // trim from left & right
  inline std::string &trim(std::string &s, const char *t = " \t\n\r\f\v")
  {
    return ltrim(rtrim(s, t), t);
  }

  std::string seekLine(std::ifstream &fs, const std::string &fieldName)
  {
    std::string line;
    std::string errMsg;
    while (std::getline(fs, line))
    {
      auto colonPos = line.find_first_of(":");
      if (colonPos == std::string::npos)
      {
        errMsg = fmt::format(
            "The TSPLIB file parser can't find [{}].",
            fieldName);
        SPDLOG_DEBUG(errMsg);
        throw std::invalid_argument(errMsg);
      }

      auto fieldNamePos = line.find(fieldName);
      if (fieldNamePos > colonPos)
        continue; // scan for next line

      // by now, we have found the field name

      // if (fieldNamePos != 0)
      // {
      //   SPDLOG_WARN(
      //     "found field name [{}] but there is preceding space(s) in that line.",
      //     fieldName);
      // }

      // grab the payload while removing all preceding and
      // trailing spaces if any (the removal is useful if T is string)
      std::string payload = line.substr(colonPos + 1);
      return trim(payload);
    }
    errMsg = fmt::format(
        "Failed to find [{}] after scanning the whole file.",
        fieldName);
    throw std::invalid_argument(errMsg);
  }

  
  Eigen::Matrix<float, -1, -1> seekSection(
      std::ifstream &fs, const std::string &sectionName,
      size_t numRows, size_t numPayloadCols, bool hasEnumerationCol)
  {
    using string = std::string;
    string line;
    string errMsg;
    if (numPayloadCols == 0)
      throw std::invalid_argument("numPayloadCols must be positive");

    while (std::getline(fs, line))
    {
      if (line.find(sectionName) == std::string::npos)
        continue;
      else // found the section, time to scan the rest
      {
        Eigen::Matrix<float, -1, -1> mat(numRows, numPayloadCols);

        for (size_t ii = 0; ii < numRows; ++ii)
        {
          if (!std::getline(fs, line))
          {
            errMsg = fmt::format(
                "Incomplete data section [{}]: fails to find row {:d} (1-based indexing)",
                sectionName, ii + 1);
            SPDLOG_ERROR(errMsg);
            throw std::invalid_argument(errMsg);
          }
          else
          {
            if (line.empty())
            {
              errMsg = fmt::format(
                  "Empty line {:d} in data section [{}]",
                  ii + 1, sectionName);
              SPDLOG_ERROR(errMsg);
              throw std::invalid_argument(errMsg);
            }
            std::stringstream ss;
            ss << line;
            std::string elem; // no space in between
            if (hasEnumerationCol)
            {
              std::getline(ss, elem, ' ');                    // a single space
              int parsedLineNumber = std::atoi(elem.c_str()); // no exception
              if (parsedLineNumber != ii + 1)
              {
                errMsg = fmt::format(
                    "Line {:d} of data section {} (see below) doesn't begin with "
                    "the expected counter:\n{}",
                    ii + 1, sectionName, line);
                SPDLOG_ERROR(errMsg);
                throw std::invalid_argument(errMsg);
              }
            }
            for (size_t jj = 0; jj < numPayloadCols; ++jj)
            {
              if (!std::getline(ss, elem, ' '))
              {
                // fewer stuff than expected
                errMsg = fmt::format(
                    "Not enough entries in line {:d} of data section {}: "
                    "expect {:d} but got only {:d}",
                    ii + 1, sectionName, numPayloadCols, jj - 1);
                SPDLOG_ERROR(errMsg);
                throw std::invalid_argument(errMsg);
              }
              else
              {
                try
                {
                  mat(ii, jj) = stof(elem); // may throw exception
                }
                catch (const std::invalid_argument &ex)
                {
                  errMsg = fmt::format(
                    "Failed to parse line {:d} because it is not entirely numbers:\n{}",
                    ii+1, line);
                  SPDLOG_ERROR(errMsg);
                  throw std::invalid_argument(errMsg);
                }
              }
            }
            if (!ss.eof())
            {
              /// more stuff there than expected
              errMsg = fmt::format(
                  "Expect only {:d} entries in line {:d} of data section {} "
                  "but got more (see below):\n{} ",
                  numPayloadCols, ii + 1, sectionName, line);
              SPDLOG_ERROR(errMsg);
              throw std::invalid_argument(errMsg);
            }
          }
          std::stringstream ss;
          ss << mat.row(ii);
          SPDLOG_DEBUG("Completed parsing row {:d} : {}", ii + 1, ss.str());
        }
        return mat;
      }
    }
    errMsg = "This part shouldn't be reached";
    SPDLOG_ERROR(errMsg);
    throw std::invalid_argument(errMsg);
  }
}