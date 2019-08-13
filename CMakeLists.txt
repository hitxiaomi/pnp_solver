cmake_minimum_required(VERSION 3.14)
project(pnp_solver)

set(DEFAULT_BUILD_TYPE "Debug")

set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")


# third party libs
# eigen
find_package(Eigen3 REQUIRED)
message(STATUS "EIGEN3_DIR:"  ${EIGEN3_INCLUDE_DIR})
include_directories(
        ${PROJECT_SOURCE_DIR}
        ${EIGEN3_INCLUDE_DIR} )

add_executable(test_dlt test_dlt.cpp pnp_solver.cpp)
add_executable(test_epnp  pnp_solver.cpp test_epnp.cpp)