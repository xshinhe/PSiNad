file(GLOB_RECURSE dev_srcs *.cpp *.c)

add_library(dev ${dev_srcs})

target_link_libraries(
    dev
    utils
    models
    base
    md
    glog::glog
    MPIMod::MPI
)

target_precompile_headers(dev REUSE_FROM utils)
