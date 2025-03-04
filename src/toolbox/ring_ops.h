#pragma once

#include <vector>
#include <cstddef>

namespace xtsp::internal
{
  /**
   * @brief Reverse the segment from [rankStart, rankEnd] (both inclusively)
   *
   * Expect rankStart < rankEnd
   */
  template <typename T>
  void reverseRingSegment_strict(std::vector<T> &ring, size_t rankStart, size_t rankEnd);

  /// @brief Automatically decide which side to flip.
  /// It's smarter than the strict version "smart" because 
  /// it can reduce the number of operations.
  /// @return true if the reversal concerns [segStart, segEnd], 
  ///    false if it's the other (geodesis) "complement" of the segment.
  template <typename T>
  bool reverseRingSegment_smart(std::vector<T> &ring, size_t segStart, size_t segEnd);
}