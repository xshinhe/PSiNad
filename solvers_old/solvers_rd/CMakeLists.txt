file(GLOB_RECURSE rd_srcs *.cpp)

add_library(rd ${rd_srcs})

target_link_libraries(
    rd
    utils
    models
    base
    md
    nad
    glog::glog
    MPIMod::MPI
)
