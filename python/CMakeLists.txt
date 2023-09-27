set(CMAKE_CXX_FLAGS "-Os -fopenmp -march=native -fPIC")

add_library(opendf SHARED opendf_test.cpp)

target_link_libraries(
    opendf
    gflags::gflags
    glog::glog
    MKLMod::MKL
    MPIMod::MPI
    HighFive
    pybind11::embed
    models
    solvers
)

