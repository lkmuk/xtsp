#include "xtsp/core/utils.h"

#include <gtest/gtest.h>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>


namespace fmt = spdlog::fmt_lib;

TEST(assertNoDuplicate, success)
{
  const std::vector<size_t> testArr {34, 23223, 123};
  xtsp::utils::assertNoDuplicate(testArr);
}

TEST(assertNoDuplicate, catchDuplicate)
{
  const std::vector<size_t> testArr {34, 23223, 123, 23223};
  const std::string arrName = "vector";
  const std::string entryName = "element";

  // change this
  const std::string expectedMsg = 
    "Invalid vector because element 23223 appears at least twice";
  try
  {
      // change this
      xtsp::utils::assertNoDuplicate(testArr, arrName, entryName);
      FAIL() << fmt::format("Should have thrown an exception: {}", expectedMsg);
  }
  catch (const std::invalid_argument& actualException)
  {
      EXPECT_EQ(actualException.what(), expectedMsg);
  }
}


TEST(assertIsPermutation, success)
{
  const std::vector<size_t> testPerm {4, 0, 1, 3, 2};
  xtsp::utils::assertIsPermutation(5, testPerm);
}

TEST(assertIsPermutation, catchMismatchSize)
{
  const std::vector<size_t> testPerm {4, 0, 1, 3, 2};
  const std::string arrName = "vector";
  const std::string entryName = "element";

  // change this
  const std::string expectedMsg = 
    "Invalid vector because of mismatched length: "
    "expect 4 got 5";
  try
  {
      // change this
      xtsp::utils::assertIsPermutation(4, testPerm, arrName, entryName);
      FAIL() << fmt::format("Should have thrown an exception: {}", expectedMsg);
  }
  catch (const std::invalid_argument& actualException)
  {
      EXPECT_EQ(actualException.what(), expectedMsg);
  }
}

TEST(assertIsPermutation, catchDuplicate)
{
  const std::vector<size_t> testPerm {4, 0, 1, 1, 2};
  const std::string arrName = "tour";
  const std::string entryName = "city";

  // change this
  const std::string expectedMsg = 
    "Invalid tour because city 1 appears at least twice "
    "(at tour positions 2 and 3.)";
  try
  {
      // change this
      xtsp::utils::assertIsPermutation(5, testPerm, arrName, entryName);
      FAIL() << fmt::format("Should have thrown an exception: {}", expectedMsg);
  }
  catch (const std::invalid_argument& actualException)
  {
      EXPECT_EQ(actualException.what(), expectedMsg);
  }
}

TEST(assertIsPermutation, catchInvalidEntry)
{
  const std::vector<size_t> testPerm {4, 0, 100, 1, 2};
  const std::string arrName = "tour";
  const std::string entryName = "city";

  // change this
  const std::string expectedMsg = 
    "Invalid tour because at position 2, "
    "the city is 100, which exceeds 4 or "
    "maybe the tour is shorter than it should";
  try
  {
      // change this
      xtsp::utils::assertIsPermutation(5, testPerm, arrName, entryName);
      FAIL() << fmt::format("Should have thrown an exception: {}", expectedMsg);
  }
  catch (const std::invalid_argument& actualException)
  {
      EXPECT_EQ(actualException.what(), expectedMsg);
  }
}