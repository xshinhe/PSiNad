cmake_minimum_required(VERSION 3.12)

file(GLOB_RECURSE ker_srcs ./*.cpp ./*.c)

add_library(ker ${ker_srcs})

set_property(TARGET ker PROPERTY POSITION_INDEPENDENT_CODE ON)

# target_precompile_headers(utils PUBLIC definitions.h phys.h)

target_link_libraries(ker
    core
    pl
)
