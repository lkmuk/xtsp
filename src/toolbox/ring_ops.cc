#include "ring_ops.h"

#include <algorithm> // std::reverse
// #include <cstring> // for memcpy
#include <stdexcept>

#include <spdlog/spdlog.h>

namespace xtsp::internal
{
  template <typename T>
  void reverseRingSegment_strict(std::vector<T> &ring, size_t rankStart, size_t rankEnd)
  {
    int segSz = (int)(rankEnd) - (int)(rankStart) + 1;
    if (segSz > (int)ring.size())
      throw std::invalid_argument("invalid segment specification");
    if (segSz <= 0)
    {
      SPDLOG_WARN("ignoring no-op request: rankStart = {:d}, rankEnd = {:d}",
        rankStart, rankEnd);
      return;
    }
    if (segSz == 1)
    {
      std::swap(ring[rankStart%ring.size()], ring[rankEnd%ring.size()]);
      return;
    }
    if (segSz == (int)ring.size())
    {
      std::reverse(std::begin(ring), std::end(ring));
      return;
    }

    // compute the "canonical" segment index
    // where the segment start index must be in {0, 1, ..., arraySz - 1}
    const size_t normalizedStart = rankStart % ring.size();
    const size_t normalizedStop = normalizedStart + segSz;
    if (normalizedStop <= ring.size())
    {
      std::reverse(
          std::begin(ring) + normalizedStart,
          std::begin(ring) + normalizedStop);
      return;
    }

    // now, we tackle the remaining case where
    // the segment wraps around the array's highest index.
    // Actually, thanks to our wrapping indexing operator,
    // we don't really need to worry about these things
    //
    // theoretical complexity of the current implementation: O(segSz)
    size_t numSwaps = (segSz >> 1); // divide by 2 --> floor

    // try setenv OMP_NUM_THREADS=4
    // #pragma omp parallel for
    for (size_t ii = 0; ii < numSwaps; ++ii)
    {
      size_t rankX = (normalizedStart + ii) % ring.size();
      size_t rankY = (normalizedStop - 1 - ii) % ring.size();
      std::swap(ring[rankX], ring[rankY]);
    }

    // Confirmed: This one isn't clearly faster (on release build)
    // for (size_t indA = normalizedStart, indB = normalizedStop - 1; indA < indB; ++indA, --indB)
    // {
    //     std::swap(this[indA%ring.size()], this[indB%ring.size()]);
    // }
  }

  template <typename T>
  bool reverseRingSegment_smart(std::vector<T> &ring, size_t segStart, size_t segEnd)
  {
    if (segStart > segEnd)
      throw std::invalid_argument("segStart index too large");
    size_t segSz = segEnd - segStart;
    if (segSz <= (ring.size() >> 1))
    {
      reverseRingSegment_strict(ring, segStart, segEnd);
      return true;
    }
    else
    {
      reverseRingSegment_strict(ring, segEnd, segEnd + ring.size() - segSz);
      return false;
    }
  }

  // template instantiation
  template bool reverseRingSegment_smart<size_t>(std::vector<size_t> &, size_t, size_t);
  template void reverseRingSegment_strict<size_t>(std::vector<size_t> &, size_t, size_t);
}
