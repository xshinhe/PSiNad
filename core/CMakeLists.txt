cmake_minimum_required(VERSION 3.12)

file(GLOB_RECURSE core_srcs ./*.cpp ./*.c)

add_library(core ${core_srcs})

set_property(TARGET core PROPERTY POSITION_INDEPENDENT_CODE ON)

# target_precompile_headers(utils PUBLIC definitions.h phys.h)

target_link_libraries(core
    fparser
    pl
)
