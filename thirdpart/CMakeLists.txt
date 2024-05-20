find_package(Git QUIET)
mark_as_advanced(CLEAR GIT_EXECUTABLE)
if(GIT_FOUND)
#  execute_process(
#    COMMAND ${GIT_EXECUTABLE} pull --recurse-submodules
#    OUTPUT_VARIABLE _dummy
#    OUTPUT_STRIP_TRAILING_WHITESPACE
#    ERROR_QUIET
#    WORKING_DIRECTORY
#      ${CMAKE_CURRENT_SOURCE_DIR}
#    )
else()
    message(STATUS "GIT is not found! Please manually download the libraries in thirdpart")
endif(GIT_FOUND)

set (BUILD_TESTING OFF)      # disable unit testing in gflags/glog/catch2
set (BUILD_SHARED_LIBS ON)   # enable gflags/glog shared library
set (WITH_GTEST OFF)

# one some old computer (old centos), after GCC5.1 it add new realization for
#     std::string & std::list in `libstdc++`. By adding `__cxx11` namespace.
#     some older compilers don't open this option by default. But it is needed
#     by gflags/glog here

add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)
add_subdirectory(gflags)
set (WITH_GFLAGS OFF)
add_subdirectory(glog)

# add_subdirectory(gperftools)
add_subdirectory(fparser)
# add_subdirectory(googletest) # depreciated; catch2 as alternative (header-only)

set (BUILD_TESTING ON)       # recover unit testing
set (BUILD_SHARED_LIBS OFF)  # disable shared library

# palanteer is a header-only profiling library.
# add_subdirectory(palanteer) # header-only!

# @depreciated
# add_definitions(-DH5_HAVE_PARALLEL=ON)
# add_definitions(-DHIGHFIVE_PARALLEL_HDF5=ON)
# add_subdirectory(HighFive)

add_subdirectory(pybind11)

