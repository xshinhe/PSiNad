set(CMAKE_CXX_FLAGS "-Os -fopenmp -march=native -fPIC")

add_library(pykids SHARED libpykids.cpp)

target_link_libraries(
    pykids
    pybind11::embed
    cppkids
)
