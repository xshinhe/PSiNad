cmake_minimum_required(VERSION 3.12)

file(GLOB_RECURSE model_srcs ./*.cpp ./*.c)

add_library(model ${model_srcs})

target_link_libraries(model
    core
    ker
    pl
)
