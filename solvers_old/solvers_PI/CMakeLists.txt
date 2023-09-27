file(GLOB_RECURSE pimd_srcs *.cpp *.c)

add_library(PI ${pimd_srcs})

target_link_libraries(
    PI
    md
    utils
    base
    models
    glog::glog
)

