#include "xtsp/core/utils.h"

#include <gtest/gtest.h>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>
#include <spdlog/fmt/bundled/ranges.h>


TEST(genPermutation, success)
{
  const size_t N = 16;
  std::vector<size_t> buf;
  xtsp::utils::Rng_T rng(123);
  xtsp::utils::genPermutation(rng, N, buf);
  ASSERT_EQ(buf.size(), N);
  xtsp::utils::assertIsPermutation(N, buf);

  spdlog::set_level(spdlog::level::info);
  SPDLOG_INFO("Generated permutation: {}",
    fmt::join(buf, " ")
  );
}