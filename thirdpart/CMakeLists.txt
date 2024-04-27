
# disable unit testing in gflags/glog/catch2
set (BUILD_TESTING OFF)
# enable gflags/glog shared library
set (BUILD_SHARED_LIBS ON)
set (WITH_GTEST OFF)
set (WITH_PKGCONFIG ON)

# one some old computer (old centos), after GCC5.1 it add new realization for
#     std::string & std::list in `libstdc++`. By adding `__cxx11` namespace.
#     some older compilers don't open this option by default. But it is needed
#     by gflags/glog here
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)
add_subdirectory(gflags)
if(NOT TARGET gflags)
    add_subdirectory(glog)
endif()
# add_subdirectory(glog)
# add_subdirectory(gperftools)
#add_subdirectory(palanteer)
#add_subdirectory(fparser)

# add_subdirectory(googletest) # depreciated; catch2 as alternative (header-only)
add_subdirectory(catch2)

# recover unit testing
set (BUILD_TESTING ON)
set (BUILD_SHARED_LIBS OFF)

# palanteer is a header-only profiling library.
# add_subdirectory(palanteer) # header-only!

# @depreciated
# add_definitions(-DH5_HAVE_PARALLEL=ON)
# add_definitions(-DHIGHFIVE_PARALLEL_HDF5=ON)
# add_subdirectory(HighFive)

add_subdirectory(pybind11)

