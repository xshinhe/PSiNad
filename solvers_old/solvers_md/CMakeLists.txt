file(GLOB_RECURSE md_srcs *.cpp *.c)

add_library(md ${md_srcs})
target_link_libraries(
    md
    utils
    base
    models
    glog::glog
    HighFive
    MKLMod::MKL
    MPIMod::MPI
)

target_precompile_headers(md REUSE_FROM utils)
