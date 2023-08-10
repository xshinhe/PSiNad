cmake_minimum_required(VERSION 3.12)

file(GLOB_RECURSE solver_srcs ./*.cpp ./*.c)

add_library(solver ${solver_srcs})

# target_precompile_headers(utils PUBLIC definitions.h phys.h)

target_link_libraries(solver
    core
    ker
    model
    pl
)
