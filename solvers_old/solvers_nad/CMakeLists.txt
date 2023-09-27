file(GLOB_RECURSE nad_srcs *.cpp *.c)

add_library(nad ${nad_srcs})

target_link_libraries(
    nad
    utils
    models
    base
    md
    glog::glog
    MPIMod::MPI
)

target_precompile_headers(nad REUSE_FROM utils)
