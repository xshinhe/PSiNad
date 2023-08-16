set(CMAKE_CXX_FLAGS "-Os -fopenmp -march=native -fPIC")

add_library(pyopendf SHARED libpyopendf.cpp)

target_link_libraries(
    pyopendf
    pybind11::embed
    models
    solvers
)
