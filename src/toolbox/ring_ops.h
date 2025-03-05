#pragma once

#include <vector>
#include <cstddef>

namespace xtsp::internal
{
  /**
   * @brief efficient inplace reversal/flip a segment indicated by [segStart, segEnd]
   * 
   * Assumption on @a segEnd and @a segStart 
   * segSz in [0, array size] where segSz =  @a segEnd - @a segStart
   * 
   * We can have both @a segStart and @a segEnd greater than the ring size
   *      
   * Special cases
   * 
   *  * with no-op, i.e., segSz = 0 or 1
   *  * segSz = the array size --- i.e., reversing the whole array
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