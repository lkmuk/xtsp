#pragma once

#include <xtsp/core/tour.h>

namespace xtsp
{
  /**
   *  Tour in adjacency table: vertex ID -> (prev vertex , next vertex) 
   */
  class AdjTabTour : public AbstractTour
  {
  public:
    /////////////////////////////
    ///  Compulsory functions
    /////////////////////////////

    /// @brief initialization from a permutation vector
    /// @param checks no revisit and no invalid index
    /// @see same API as PermTour::PermTour
    AdjTabTour(const std::vector<size_t>& perm, int maxSize = -1, bool checks = true);

    virtual size_t size() const override;
    virtual size_t maxSize() const override;
    virtual size_t next(size_t vertex) const override;
    virtual size_t getDepotId() const override;

    virtual void exchangeTwoEdges(
      size_t vA, size_t vC, bool strict = false) override;

    /// @todo add to Public API
    // virtual bool isBetween(size_t vA, size_t vB, size_t vC) const override;

    //////////////////////////////////////////////
    /// Specialized implementation of public API
    /////////////////////////////////////////////
    virtual bool isOneStepAhead(size_t vStart, size_t vGoal) const override
    {
      return m_dat[vStart].next == &m_dat[vGoal];
    }
    virtual bool isTwoPlusStepsAhead(size_t vStart, size_t vGoal) const override
    {
      return !((vStart==vGoal)||isOneStepAhead(vStart, vGoal));
    }

    ///////////////////////////////
    ///  Non-standard API functions
    ///////////////////////////////

    // size_t prev(size_t vertex) const
    // {
    //   return m_dat[vertex].prev->id;
    // }
    
    // // manipulate the tour
    // void setDepot(size_t vDepot)
    // {
    //   m_head = vDepot;
    // }
    // void reverse() 
    // {
    //   for (auto& n : m_dat)
    //   {
    //     std::swap(n.prev, n.next);
    //   }
    // }

    // 
    
    
  protected:
    struct Vertex
    {
      const size_t id;  // ID of this vertex
      Vertex* prev = nullptr;
      Vertex* next = nullptr;
      // In our case, the ID of this vertex is already 
      // implied by the address (&this) but it is included 
      // anyways to simplify the implementation.

      Vertex(size_t v) : id(v)
      {
      }
    };
    // the adjacency table (which shouldn't get expanded)
    std::vector<Vertex> m_dat; 
    size_t m_head = 0; // used for iterator
    size_t m_cache_tourSize = 0;
    const size_t m_N; // maximum size
  };

}