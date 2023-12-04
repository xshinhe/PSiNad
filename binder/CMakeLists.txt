set(CMAKE_CXX_FLAGS "-Os -fopenmp -march=native -fPIC")

add_library(opendfpy SHARED libopendfpy.cpp)

target_link_libraries(
    opendfpy
    pybind11::embed
    opendfcpp
)
