#ifndef __XTSP_CORE_UTILS__
#define __XTSP_CORE_UTILS__

#include <vector>
#include <cstddef>
#include <string>
#include <random>

namespace xtsp
{
  namespace utils
  {
    // We commit to this RNG algorithm
    using Rng_T = std::mt19937;
    
    /// generate a random permutation from 0 to N-1
    void genPermutation(Rng_T rng, size_t N, std::vector<size_t>& out);

    /// ensure every entry of an array is unique
    /// @param arr
    /// @param[in] arrName how to call an array 
    ///        in the error message, e.g., "array".
    /// @param[in] entryName how to call an entry of a permutation 
    ///        in the error message, e.g., "element".
    void assertNoDuplicate(
      const std::vector<size_t> vec,
      const std::string arrName = "array",
      const std::string entryName = "element");

    /// @brief assert a permutation vector is valid
    /// @param[in] N is the expected length of the permutation (which can be zero)
    /// @param[in] perm a permuation from 0 up to N-1
    /// @param[in] permName how to call a permutation 
    ///        in the error message, e.g., "tour".
    /// @param[in] entryName how to call an entry of a permutation 
    ///        in the error message, e.g., "city".
    void assertIsPermutation(
      size_t N,
      const std::vector<size_t> perm,
      const std::string permName = "permutation", // or "tour"
      const std::string entryName = "vertex" // or "city"
      );

  } // namespace utils

} // namespace xtsp

#endif