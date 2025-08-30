cmake_minimum_required(VERSION 3.17)
project(test_kernel_naforce LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use the library and includes from the main project
include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/core/include
  ${CMAKE_SOURCE_DIR}/src/kernels/include
  ${CMAKE_SOURCE_DIR}/src/models/include
  ${CMAKE_SOURCE_DIR}/src/solvers/include
  ${CMAKE_SOURCE_DIR}/thirdpart
  ${CMAKE_SOURCE_DIR}/thirdpart/Eigen
)

add_executable(test_kernel_naforce main.cpp)

# Link against the shared library built by the main project
# Assumes you configured and built the top-level with PSINAD_BUILD_SHARED_LIB=ON
find_library(PSINAD_LIB NAMES PSiNad_shared PATHS ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/.. ${CMAKE_SOURCE_DIR}/build)
if(NOT PSINAD_LIB)
  message(FATAL_ERROR "Could not find PSiNad_shared library. Build the top-level project first.")
endif()

# Also link gflags and glog built by the main project (find them by names as they are built locally)
find_library(GFLAGS_LIB NAMES gflags PATHS ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/build)
find_library(GLOG_LIB   NAMES glog   PATHS ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/build)

if(GFLAGS_LIB)
  target_link_libraries(test_kernel_naforce PRIVATE ${GFLAGS_LIB})
endif()
if(GLOG_LIB)
  target_link_libraries(test_kernel_naforce PRIVATE ${GLOG_LIB})
endif()

target_link_libraries(test_kernel_naforce PRIVATE ${PSINAD_LIB})
