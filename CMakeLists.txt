cmake_minimum_required(VERSION 3.12)
project(xtsp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_POSITION_INDEPENDENT_CODE FALSE)

# message(${CMAKE_CXX_COMPILER_ID})
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
message("setting extra compiler flags")    
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -ffast-math")
endif()

find_package(Eigen3 CONFIG REQUIRED)
find_package(GTest CONFIG REQUIRED)
add_subdirectory(deps/spdlog)

add_library(${PROJECT_NAME}
    src/core/utils.cc
    src/core/clustering.cc
    src/core/complete_graph.cc
    src/core/tour.cc
    src/core/tsplib_io.cc
    src/core/tsplib_io_seek_impl.cc
    src/initialization/insertion.cc
    src/local_search/kopt.cc
    src/local_search/gtsp_only.cc
    src/toolbox/ring_ops.cc
)
target_include_directories(${PROJECT_NAME}
    PUBLIC include 
)
target_link_libraries(${PROJECT_NAME}
    PUBLIC Eigen3::Eigen
    PUBLIC spdlog::spdlog
)

add_executable(demo_load_geom_tsp
  demo/demo_load_geom_TSP.cc
)
target_link_libraries(demo_load_geom_tsp ${PROJECT_NAME})

enable_testing()

add_executable(test_core
    tests/core/test_utils.cc
    tests/core/test_utils_rng.cc
    tests/core/test_clustering.cc
    tests/core/test_complete_graph.cc
    tests/core/test_tour.cc
)
target_link_libraries(test_core PRIVATE ${PROJECT_NAME} GTest::gtest GTest::gtest_main)
add_test(test_core ${CMAKE_BINARY_DIR}/test_core)

# better to split out these test cases 
# so that we can focus on them if we need to.
add_executable(test_io
    tests/core/test_tsplib_io.cc
)
target_link_libraries(test_io PRIVATE ${PROJECT_NAME} GTest::gtest GTest::gtest_main)
add_test(test_io ${CMAKE_BINARY_DIR}/test_io)

add_executable(test_constructTour
    tests/initialization/test_insertion.cc
)
target_link_libraries(test_constructTour 
  PRIVATE ${PROJECT_NAME} GTest::gtest GTest::gtest_main)
add_test(test_constructTour ${CMAKE_BINARY_DIR}/test_constructTour)

add_executable(test_local_search
    tests/local_search/test_gtsp_only.cc
    tests/local_search/test_kopt.cc
)
target_link_libraries(test_local_search PRIVATE ${PROJECT_NAME} GTest::gtest GTest::gtest_main)
add_test(test_local_search ${CMAKE_BINARY_DIR}/test_local_search)
