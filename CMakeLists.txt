cmake_minimum_required(VERSION 3.16)

# choose a different directory from source for compilation
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "
    FATAL: In-source builds are not allowed.
           You should create a separate directory for build files.
           ")
endif()

project(opendf LANGUAGES C CXX Fortran)
enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/thirdpart/catch2/contrib)

# processing repositary information
include(repo)
get_repo_info(GIT_BRANCH GIT_DATE GIT_HASH GIT_STATUS)
configure_file(
  ${CMAKE_SOURCE_DIR}/version/version.h.in
  ${CMAKE_SOURCE_DIR}/generate/version.h
  @ONLY
)

# processing modular configuration (solver / model)
# include(module)
# module_config()

# set ccahe to accelerate compile
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    message(STATUS "set up ccache ...")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif()

include_directories(${CMAKE_SOURCE_DIR}/thirdpart/)
include_directories(${CMAKE_SOURCE_DIR}) # set up root also as include directory

################################################################################
# set esssential variables
################################################################################


option (BUILD_PYTHON_API "Build Python API" OFF)
option (BUILD_SHARED_LIBS "Build shared libraries" OFF)

# set(CMAKE_BUILD_TYPE <Config>)
## Please -DCMAKE_BUILD_TYPE=<Config> : Debug, Release, RelWithDebInfo, MinsizeRel
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -march=native -fPIC")

# ################################################################################
# # configure libraries
# ################################################################################
# set(CUSTOMIZED_MKL_DIR "/opt/intel/oneapi/mkl/latest/")
# set(CUSTOMIZED_MPI_DIR "/opt/intel/oneapi/mpi/latest/")
# set(CUSTOMIZED_HDF5_DIR "/usr/local/HDF_Group/HDF5/1.13.2")

find_package(MPIMod REQUIRED)
# find_package(MKLMod)

# find_package(CINTMod) # optional
# # find_package(pybind11 REQUIRED)

# #set(AMBERHOME /home/zhf/work/amberfold/amber20/) # @blame DON'T WRITE YOUR LOCAL DIR HERE!

# include_directories(${AMBERHOME}/include/)
# link_directories(${AMBERHOME}/lib/)

# ################################################################################
# # compilation of different parts
# ################################################################################

# # thirdpart library
add_subdirectory(thirdpart)

add_subdirectory(core)
add_subdirectory(kernels)
add_subdirectory(models)
add_subdirectory(solvers)
# add_subdirectory(python)

add_executable(opendf main.cc Handler.cpp mpi_utils.cpp)
target_link_libraries(opendf
    # models
    # solvers
    gflags::gflags
    glog::glog
    fparser
    # MKLMod::MKL
    MPIMod::MPI
    core
    ker
    model 
    solver
    pl
)

set(OPENDF_UNIT_TESTS ON)
if(OPENDF_UNIT_TESTS)
  enable_testing()
  # add_subdirectory(tests)
endif()
