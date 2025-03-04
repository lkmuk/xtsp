#pragma once

#include <vector>
#include <cstddef>

#include <algorithm>
#include <spdlog/spdlog.h>

namespace xtsp::internal
{
  /**
   * Keep track of an On/Off state of N items
   *
   * We will identify each item with a unsigned
   * integer number.
   *
   *  Operations
   *
   *  A. query if there is still elements to be processed
   *  B. random access
   *  C. remove an element from the buffer
   *
   *  *  Possible data structures and its worst-case complexity
   *
   *  1. std::set
   *     (A: O(log N), B: O(log N), C: O(log N))
   *  2. bitfield using std::vector<bool>
   *     (A: O(N), B : O(1), C: O(1))
   *
   * Theoretically, std::set has more scalable query operation.
   * However, typical usage suggest bitfields are often superior
   * in runtime.
   *
   * Despite that, using std::vector<bool> directly can be a bit
   * unintuitive or even error-prone, hence this helper class.
   */
  class WorkBuffer
  {
  public:
    // initialize all items as TODO
    WorkBuffer(size_t numItems) : m_dat(numItems, true)
    {
    }
    bool isEmpty() const
    {
      return (std::find(m_dat.cbegin(), m_dat.cend(), true) == m_dat.cend());
    }
    bool isDone(size_t queryItem) const
    {
      return m_dat[queryItem];
    }
    // not really a tight bound ...
    size_t estimateTodoIdMin() const
    {
      return 0;
    }
    // not really a tight bound (inclusive)
    size_t estimateTodoIdMax() const
    {
      return m_dat.size()-1;
    }


    void markDone(size_t item)
    {
      m_dat[item] = false;
    }
    void markTodo(size_t item)
    {
      m_dat[item] = true;
    }

  protected:
    std::vector<size_t> m_dat;
  };

  // more efficient query than WorkBuffer
  // now the query operation becomes O(1)
  // but does it make sense at all?
  // 
  // result: not much as far as I can tell
  /*
  class CachedWorkBuffer
  {
  public:
    // initialize all items as TODO
    CachedWorkBuffer(size_t numItems) : m_dat(numItems, true), m_cache_max(numItems - 1)
    {
    }
    bool isEmpty() const
    {
      return m_cache_min > m_cache_max;
    }
    bool isDone(size_t queryItem) const
    {
      return m_dat[queryItem];
    }
    // the exact one (inclusive)
    size_t estimateTodoIdMin() const
    {
      return m_cache_min;
    }
    // the exact one (inclusive)
    size_t estimateTodoIdMax() const
    {
      return m_cache_max;
    }

    void markDone(size_t item)
    {
      if (m_dat[item]) // a new update
      {
        m_dat[item] = false;

        if ((int)item == m_cache_min)
        {
          m_cache_min = std::distance(
              m_dat.cbegin(),
              std::find(m_dat.cbegin(), m_dat.cend(), true));
        }
        if ((int)item == m_cache_max)
        {
          // this one is a bit technique and less intuitive
          // notice the -1 (because here, max is inclusive)
          m_cache_max = (int)std::distance(
                            std::find(m_dat.crbegin(), m_dat.crend(), true),
                            m_dat.crend()) - 1;
        }

        SPDLOG_DEBUG(
            "dequeuing {:d}: cached TODO IDs ranges from [{:d},{:d}]",
            item, m_cache_min, m_cache_max);
      }
    }
    void markTodo(size_t item)
    {
      if (!m_dat[item]) // a new update
      {
        m_dat[item] = true;
        m_cache_max = (m_cache_max < (int)item) ? item : m_cache_max;
        m_cache_min = ((int)item < m_cache_min) ? item : m_cache_min;
      }
    }

  protected:
    std::vector<size_t> m_dat; // True = TODO, false = DONE
    int m_cache_min = 0;    // min ID that is a TODO
    int m_cache_max;        // the max ID that is a TODO
  };
  */

} // namespace xtsp::internal