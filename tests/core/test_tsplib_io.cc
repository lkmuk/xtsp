// TODO: directly test the parser
#include "xtsp/core/tsplib_io.h"
// for now test this wrapper first
#include "xtsp/core/complete_graph.h"
// we also want to test tour file I/O
#include "xtsp/core/tour.h" 

#include <gtest/gtest.h>
#include <filesystem>

// just in case we want to adjust the logging level
#include <spdlog/spdlog.h>

static const auto dataDir = std::filesystem::path(
  __FILE__).parent_path().parent_path()/"dataset";

TEST(tsplibParser, euc2dTspPr144)
{
  // this one actually has XY in integer
  const auto testFile = dataDir/"pr144.tsp";
  
  spdlog::set_level(spdlog::level::warn);

  auto g = xtsp::ImplicitCompleteGraph<float>::loadFromTsplibFile(testFile.c_str());

  ASSERT_EQ(g.numVertices(), 144);
  ASSERT_EQ(g.nDim(), 2);
  const Eigen::Matrix<float, -1, 2>& xyView = g.getXy();
  EXPECT_FLOAT_EQ(xyView(0,0), 4350);  EXPECT_FLOAT_EQ(xyView(0,1), 4425);
  EXPECT_FLOAT_EQ(xyView(2,0), 4300);  EXPECT_FLOAT_EQ(xyView(2,1), 4725);
  EXPECT_FLOAT_EQ(xyView(143,0), 15225); EXPECT_FLOAT_EQ(xyView(143,1), 3150);

  EXPECT_FALSE(g.isClustered());

}

TEST(tsplibParser, euc2dTspU1817)
{
  // this one actually has XY in Floating point (scientific notation)
  const auto testFile = dataDir/"u1817.tsp";
  
  spdlog::set_level(spdlog::level::warn);

  auto g = xtsp::ImplicitCompleteGraph<float>::loadFromTsplibFile(testFile.c_str());

  ASSERT_EQ(g.numVertices(), 1817);
  ASSERT_EQ(g.nDim(), 2);
  const Eigen::Matrix<float, -1, 2>& xyView = g.getXy();
  EXPECT_FLOAT_EQ(xyView(0,0), 6.5119e2);  EXPECT_FLOAT_EQ(xyView(0,1), 2.24439e3);
  EXPECT_FLOAT_EQ(xyView(1,0), 6.766e2);  EXPECT_FLOAT_EQ(xyView(1,1), 2.1682e3);
  EXPECT_FLOAT_EQ(xyView(1816,0), 6.766e2); EXPECT_FLOAT_EQ(xyView(1816,1), 2.24439e3);

  EXPECT_FALSE(g.isClustered());

}

TEST(tsplibParser, euc2dGeneralizedTsp20kroA100)
{

  const auto testFile = dataDir/"20kroA100.gtsp";

  spdlog::set_level(spdlog::level::warn);

  auto g = xtsp::ImplicitCompleteGraph<float>::loadFromTsplibFile(testFile.c_str());

  ASSERT_EQ(g.numVertices(), 100);
  ASSERT_EQ(g.nDim(), 2);
  const Eigen::Matrix<float, -1, 2>& xyView = g.getXy();
  EXPECT_FLOAT_EQ(xyView(0,0), 1380);  EXPECT_FLOAT_EQ(xyView(0,1), 939);
  EXPECT_FLOAT_EQ(xyView(99,0), 3950); EXPECT_FLOAT_EQ(xyView(99,1), 1558);

  ASSERT_TRUE(g.isClustered());
  ASSERT_EQ(g.numClusters(), 20);
  const auto clustering = g.getClusteringInfo();
  EXPECT_EQ(clustering->evalWhichHasTheLeastVertices(), 18);
  // check cluster 0 (btw everything is 0-indexed)
  ASSERT_EQ(clustering->getClusterSize(0), 4);
  EXPECT_EQ(clustering->getMembers(0)[0], 40); 
  EXPECT_EQ(clustering->getMembers(0)[1], 47);
  EXPECT_EQ(clustering->getMembers(0)[2], 70);
  EXPECT_EQ(clustering->getMembers(0)[3], 99);
  // check cluster M-1 = 19
  ASSERT_EQ(clustering->getClusterSize(19), 3);
  EXPECT_EQ(clustering->getMembers(19)[0], 29); 
  EXPECT_EQ(clustering->getMembers(19)[1], 38);
  EXPECT_EQ(clustering->getMembers(19)[2], 95);
}


/**
 * robustness and error handling 
 * 
 * There are too many possible variations (e.g., padding spaces around colons).
 * So are the kinds of errors.
 * Many were once tested but not yet implemented as unit tests yet.
 * 
 *   * file not found
 *   * can't find certain lines (e.g. NAME, TYPE, ...)
 *   * found certain lines, e.g., TYPE, EDGE_WEIGHT_TYPE but they are not supported.
 *   * can't find certain section (e.g., NODE_COORD_SECTION)
 *   * found NODE_COORD_SECTION but failed because ...
 * 
 *     * too short
 *     * failed to parse a line
 *       * cannot be converted to a number (e.g., "2 c34 4,1")
 *         (the parser will try to recover e.g., "3 4x1 5.2\n" as "3. 4. 5.2")
 *       * wrong line number (e.g., "32 12.2 -44.1" for supposedly line 3)
 *       * too many/few entries (e.g., for EUC3D problem, we expect 3 payload entries)
 *     * (too long) <-- this one we will warn the user without raising an exception
 *   
 *   * is GTSP but can't found
 *   * found GTSP_SET_SECITON but 
 *     * ...
 */
TEST(tsplibParser, catchBadCoordinateLineNumber)
{
  const auto testFile = dataDir/"malformed"/"pr10_wrongDataSectionNumber.tsp";

  spdlog::set_level(spdlog::level::info);

  const std::string expectedMsg = 
    "Line 3 of data section NODE_COORD_SECTION (see below) "
    "doesn't begin with the expected counter:\n"
    "4 4300 4725";
    
  try
  {
      // change this
      auto g = xtsp::ImplicitCompleteGraph<float>::loadFromTsplibFile(testFile.c_str());
      FAIL() << fmt::format("Should have thrown an exception: {}", expectedMsg);
  } 
  catch (const std::invalid_argument& actualException)
  {
      EXPECT_EQ(actualException.what(), expectedMsg);
  }

}

TEST(tsplibParser, catchBadCoordinateLineExtraEntry)
{
  const auto testFile = dataDir/"malformed"/"badEuc2d_badLineWithTooManyEntries.tsp";

  spdlog::set_level(spdlog::level::info);

  const std::string expectedMsg = 
    "Expect only 2 entries in line 3 of data section NODE_COORD_SECTION "
    "but got more (see below):\n"
    "3 4300 4825 32 ";
  try
  {
      // change this
      auto g = xtsp::ImplicitCompleteGraph<float>::loadFromTsplibFile(testFile.c_str());
      FAIL() << fmt::format("Should have thrown an exception: {}", expectedMsg);
  } 
  catch (const std::invalid_argument& actualException)
  {
      EXPECT_EQ(actualException.what(), expectedMsg);
  }

}

class TsplibWriteTour : public testing::Test
{
protected:
  void SetUp() override
  {
    spdlog::set_level(spdlog::level::info);
    m_tmpOutputPath = std::filesystem::temp_directory_path()/"dummyTour.tour";
    SPDLOG_INFO("tour output filepath : {}", m_tmpOutputPath.c_str());
  }

  std::filesystem::path m_tmpOutputPath;
  // this is actually a generalized tour
  std::vector<size_t> m_perm = {23, 0, 4, 2, 5}; 
};

TEST_F(TsplibWriteTour, permTour)
{
  std::string expectedString = 
    "NAME : A dummy generalized tour used for unit testing\n"
    "TYPE : TOUR\n"
    "DIMENSION : 5\n"
    "TOUR_SECTION\n"
    "24\n"
    "1\n"
    "5\n"
    "3\n"
    "6\n"
    "-1\n"
    "EOF\n";

  auto tour = xtsp::PermTour(m_perm, 50); // 50 is just an arbitrary number
  // Subject under test
  tour.saveTsplib(m_tmpOutputPath.c_str(), "A dummy generalized tour used for unit testing");
  std::ifstream retrievedStream (m_tmpOutputPath.c_str());

  std::ostringstream sstr;
  sstr << retrievedStream.rdbuf();
  std::string retrievedString = sstr.str();

  EXPECT_EQ(retrievedString, expectedString);

}