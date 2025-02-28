cmake_minimum_required(VERSION 3.12)
project(xtsp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_POSITION_INDEPENDENT_CODE FALSE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -fopenmp")

find_package(Eigen3 CONFIG REQUIRED)
find_package(GTest CONFIG REQUIRED)
add_subdirectory(deps/spdlog)

add_library(${PROJECT_NAME}
    src/core/clustering.cc
    src/core/complete_graph.cc
)
target_include_directories(${PROJECT_NAME}
    PUBLIC include 
)
target_link_libraries(${PROJECT_NAME}
    PUBLIC Eigen3::Eigen
    PUBLIC spdlog::spdlog
)


enable_testing()

add_executable(test_core
    tests/core/test_clustering.cc
    tests/core/test_complete_graph.cc
)
target_link_libraries(test_core PRIVATE ${PROJECT_NAME} GTest::gtest GTest::gtest_main)
add_test(test_core ${CMAKE_BINARY_DIR}/test_core)
