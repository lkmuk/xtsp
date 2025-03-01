#include "xtsp/core/utils.h"

#include <limits>
#include <exception>
#include <set>

#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>

namespace xtsp
{
  namespace utils
  {
    namespace fmt = spdlog::fmt_lib;

    void assertNoDuplicate(
      const std::vector<size_t> vec,
      const std::string arrName,
      const std::string entryName)
    {
      std::set<size_t> elemSet;
      for (const size_t val : vec)
      {
        const auto &[ite, isDistinct] = elemSet.emplace(val);
        if (!isDistinct)
        {
          std::string errMsg = fmt::format(
            "Invalid {} because {} {:d} appears at least twice",
            arrName, entryName, val);
          throw std::invalid_argument(errMsg);
        }
      }

      // if (elemSet.size() != vec.size())
      // {
      //   std::string errMsg = fmt::format(
      //       "Found {:d} duplicates", vec.size() - elemSet.size());
      //   throw std::invalid_argument(errMsg);
      // }
    }

    void assertIsPermutation(
        size_t N,
        const std::vector<size_t> perm,
        const std::string permName,
        const std::string entryName)
    {
      std::string errMsg;

      if (N != perm.size())
      {
        errMsg = fmt::format(
          "Invalid {} because of mismatched length: "
          "expect {:d} got {:d}",
          permName, N, perm.size());
        throw std::invalid_argument(errMsg);
      }
      const size_t uninitializedVal = std::numeric_limits<size_t>::max();
      /// Reverse Lookup: revLut = argsort(perm).
      /// i.e., revLut[i] = position of i in \p perm .
      std::vector<size_t> revLut(N, uninitializedVal);
      for (size_t rank = 0; rank < N; ++rank)
      {
        size_t thisVertex = perm[rank];
        if (thisVertex >= N)
        {
          std::string errMsg = spdlog::fmt_lib::format(
              "Invalid {} because at position {:d}, "
              "the {} is {:d}, which exceeds {:d} or "
              "maybe the {} is shorter than it should",
              permName, rank, entryName, thisVertex, N - 1, permName);
          throw std::invalid_argument(errMsg);
        }

        if (revLut[thisVertex] != uninitializedVal)
        {
          std::string errMsg = fmt::format(
              "Invalid {} because {} {:d} appears at least twice "
              "(at {} positions {:d} and {:d}.)",
              permName, entryName, thisVertex,
              permName, revLut[thisVertex], rank);
          throw std::invalid_argument(errMsg);
        }
        revLut[thisVertex] = rank;
      }

      // ensure no vertex is un-assigned
      /// @todo: remove this because this is probably not needed
      // auto iteVertex = std::find(revLut.cbegin(), revLut.cend(), uninitializedVal);
      // if (iteVertex != perm.cend())
      // {
      //   std::string errMsg = spdlog::fmt_lib::format(
      //       "Invalid {} because {} {:d} is not in the {}.",
      //       permName, entryName,
      //       std::distance(revLut.cbegin(), iteVertex),
      //       permName);
      //   throw std::invalid_argument(errMsg);
      // }
    }
  } // namespace utils
} // namespace xtsp
