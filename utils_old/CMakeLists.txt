cmake_minimum_required(VERSION 3.12)
project(utils)

file(GLOB_RECURSE utils_srcs ./*.cpp ./*.c)

add_library(utils ${utils_srcs})

target_link_libraries(
    utils
    gflags::gflags
    glog::glog
)

target_precompile_headers(utils PUBLIC definitions.h phys.h)
