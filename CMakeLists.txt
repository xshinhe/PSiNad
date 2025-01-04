cmake_minimum_required(VERSION 3.17)

cmake_policy(SET CMP0003 NEW)  # Libraries linked via full path no longer produce linker search paths.
cmake_policy(SET CMP0005 NEW)  # Preprocessor definition values are now escaped automatically.
cmake_policy(SET CMP0011 NEW)  # Included scripts do automatic cmake_policy PUSH and POP.

# Define library name & version
set(PSINAD_LIBRARY_NAME PSiNad)
if(NOT PROJECT_NAME)
    project(${PSINAD_LIBRARY_NAME} LANGUAGES C CXX)
endif()
set(PSINAD_MAJOR_VERSION 0)
set(PSINAD_MINOR_VERSION 1)
set(PSINAD_BUILD_VERSION 1)
set(PSINAD_VERSION "${PSINAD_MAJOR_VERSION}.${PSINAD_MINOR_VERSION}.${PSINAD_BUILD_VERSION}"
     CACHE STRING "This is the version of PSiNad which will be built." FORCE)
mark_as_advanced(PSINAD_VERSION)
add_definitions(-DPSINAD_LIBRARY_NAME=${PSINAD_LIBRARY_NAME}
                -DPSINAD_MAJOR_VERSION=${PSINAD_MAJOR_VERSION}
                -DPSINAD_MINOR_VERSION=${PSINAD_MINOR_VERSION}
                -DPSINAD_BUILD_VERSION=${PSINAD_BUILD_VERSION})

# Set the C++ standard to 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Python interpreter
find_package(PythonInterp REQUIRED)
mark_as_advanced(CLEAR PYTHON_EXECUTABLE)

# Prevent in-source builds
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "In-source builds are not allowed. Please use a separate build directory.")
endif()

# Option to enable verbose
option(ENABLE_VERBOSE_MAKEFILE "Enable verbose output during build" OFF)
if(ENABLE_VERBOSE_MAKEFILE)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

# Option to enable ccache
option(CCACHE_ENABLE "Enable ccache for faster compilation" OFF)
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_ENABLE AND CCACHE_PROGRAM)
    message(STATUS "ccache is enabled and found at: ${CCACHE_PROGRAM}. Setting up ccache for compilation and linking...")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ${CCACHE_PROGRAM})
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ${CCACHE_PROGRAM})
elseif(CCACHE_ENABLE)
    message(WARNING "ccache is enabled but not found. Skipping ccache setup.")
else()
    message(STATUS "ccache is disabled. Skipping ccache setup.")
endif()

# Find and enable compile-time-perf if available
find_package(compile-time-perf QUIET)
if(compile-time-perf_FOUND)
    message(STATUS "compile-time-perf found. Enabling compile-time performance measurements...")
    enable_compile_time_perf(kids-ctp)
else()
    message(STATUS "compile-time-perf not found. Skipping compile-time performance measurements.")
endif()

# Append additional module paths & include some customized modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/thirdpart/catch2/contrib")
include(printings)
include(TargetArch)
# include(module)
# module_config()

# Options to determine what type of libraries to build
option(PSINAD_BUILD_SHARED_LIB "Build shared PSiNad libraries" ON)
option(PSINAD_BUILD_STATIC_LIB "Build static PSiNad libraries" ON)

# Add include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/generate
    ${CMAKE_SOURCE_DIR}/thirdpart
    ${CMAKE_SOURCE_DIR}/thirdpart/Eigen
    ${CMAKE_SOURCE_DIR}/thirdpart/filesystem/include
)

################################################################################
# set esssential variables
################################################################################

# set(CMAKE_BUILD_TYPE <Config>)
## Please -DCMAKE_BUILD_TYPE=<Config> : Debug, Release, RelWithDebInfo, MinsizeRel
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -march=native")

if(APPLE) # disable some warnings
    add_compile_options(-Wno-switch -Wno-braced-scalar-init -Wno-\#pragma-messages -Wno-enum-compare-switch -Wno-uninitialized)
endif()

# ################################################################################
# # configure libraries
# ################################################################################
# set(CUSTOMIZED_MKL_DIR "/opt/intel/oneapi/mkl/latest/")
# set(CUSTOMIZED_MPI_DIR "/opt/intel/oneapi/mpi/latest/")
# set(CUSTOMIZED_HDF5_DIR "/usr/local/HDF_Group/HDF5/1.13.2")

#find_package(MKLMod)

# find_package(CINTMod) # optional

# #set(AMBERHOME /home/zhf/work/amberfold/amber20/) # @blame DON'T WRITE YOUR LOCAL DIR HERE!
# include_directories(${AMBERHOME}/include/)
# link_directories(${AMBERHOME}/lib/)

# ################################################################################
# # compilation of different parts
# ################################################################################

# Set default installation prefix if not provided
if(NOT DEFINED PSINAD_INSTALL_PREFIX)
    if(WIN32)
        set(PSINAD_INSTALL_PREFIX "$ENV{ProgramFiles}/PSINAD")
    else()
        set(PSINAD_INSTALL_PREFIX "/usr/local/psinad")
    endif()
endif()

# Set CMAKE_INSTALL_PREFIX if it has not been initialized by the user
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${PSINAD_INSTALL_PREFIX} CACHE PATH "Where to install PSiNad" FORCE)
endif()
# Set default build type to Release for UNIX systems if not specified
if(UNIX AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build: Debug or Release" FORCE)
endif()
# Set C++ compiler flags for Debug build type if not already set
if(NOT CMAKE_CXX_FLAGS_DEBUG)
    set(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "Flags for Debug build" FORCE)
endif()
# Set C++ compiler flags for Release build type if not already set
if(NOT CMAKE_CXX_FLAGS_RELEASE)
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "Flags for Release build" FORCE)
endif()

# Collect up information about the version of the PSINAD library we're building
# and make it available to the code so it can be built into the binaries.

# Processing repository fingerprint
include(repo)
get_repo_info(GIT_BRANCH GIT_DATE GIT_HASH GIT_STATUS)
configure_file(
  ${CMAKE_SOURCE_DIR}/version/version.h.in
  ${CMAKE_SOURCE_DIR}/src/core/include/kids/version.h
  @ONLY
)
set(PSINAD_GIT_VERSION ${GIT_HASH})

# Define library targets
#set(SHARED_TARGET ${PSINAD_LIBRARY_NAME}) # revise to the following line
set(SHARED_TARGET ${PSINAD_LIBRARY_NAME}_shared)
set(STATIC_TARGET ${PSINAD_LIBRARY_NAME}_static)

# If no custom executable output path is specified, set it to the binary directory
if(NOT EXECUTABLE_OUTPUT_PATH)
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}
        CACHE INTERNAL "Single output directory for building all executables.")
endif()

# If no custom library output path is specified, set it to the binary directory
if(NOT LIBRARY_OUTPUT_PATH)
    set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}
        CACHE INTERNAL "Single output directory for building all libraries.")
endif()

# Define variables for executable and library directories
set(${PROJECT_NAME}_EXECUTABLE_DIR ${EXECUTABLE_OUTPUT_PATH}/${CMAKE_CFG_INTDIR})
set(${PROJECT_NAME}_LIBRARY_DIR    ${LIBRARY_OUTPUT_PATH}/${CMAKE_CFG_INTDIR})


# The source is organized into subdirectories, but we handle them all together
# rather than visit them as SUBDIRS.
set(PSINAD_SOURCE_SUBDIRS
    src/core
    src/kernels
    src/models
    src/solvers
)
file(GLOB third_src_files  ${CMAKE_SOURCE_DIR}/thirdpart/fparser/*.cpp)

# Set up API include directories
set(API_INCLUDE_DIRS) # start empty
foreach(subdir ${PSINAD_SOURCE_SUBDIRS})
    list(APPEND API_INCLUDE_DIRS
        ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include/kids
    )
endforeach()

# Collect absolute and relative include files
set(API_REL_INCLUDE_FILES)
set(API_ABS_INCLUDE_FILES)
foreach(dir ${API_INCLUDE_DIRS})
    file(GLOB_RECURSE fullpaths ${dir}/*.h) # returns full pathnames
    list(APPEND API_ABS_INCLUDE_FILES ${fullpaths})
    foreach(pathname ${fullpaths})
        get_filename_component(filename ${pathname} NAME)
        list(APPEND API_REL_INCLUDE_FILES ${dir}/${filename})
    endforeach()
endforeach()

# Collect source files
set(SOURCE_FILES)
set(SOURCE_INCLUDE_FILES)
foreach(subdir ${PSINAD_SOURCE_SUBDIRS})
    file(GLOB_RECURSE src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp)
    file(GLOB_RECURSE incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    list(APPEND SOURCE_FILES ${src_files})
    list(APPEND SOURCE_INCLUDE_FILES ${incl_files})
    include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
endforeach()

# Set compile flags for specific source files
set_source_files_properties(
    ${CMAKE_SOURCE_DIR}/src/core/src/linalg.cpp
    ${CMAKE_SOURCE_DIR}/src/core/src/linalg_1.cpp
    PROPERTIES COMPILE_FLAGS "-DHAVE_SSE2=1"
)

# Include directories before compiling source files
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(EXTRA_LINK_FLAGS ${EXTRA_COMPILE_FLAGS})
if (NOT MSVC)
    if (CMAKE_SYSTEM_NAME MATCHES "Linux" AND NOT ANDROID)
        set(EXTRA_LINK_FLAGS "${EXTRA_LINK_FLAGS} -Wl,--no-as-needed -lrt")
    endif()
endif()

set(SOURCE_FILES ${SOURCE_FILES} ${third_src_files})

if(PSINAD_BUILD_SHARED_LIB)
    add_library(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})
    set_target_properties(${SHARED_TARGET} PROPERTIES 
        LINK_FLAGS "${EXTRA_LINK_FLAGS}" 
        COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DPSINAD_BUILDING_SHARED_LIBRARY -DPTHREAD_BUILDING_SHARED_LIBRARY" 
        SOVERSION "${PSINAD_MAJOR_VERSION}.${PSINAD_MINOR_VERSION}")
endif()

if(PSINAD_BUILD_STATIC_LIB)
    add_library(${STATIC_TARGET} STATIC ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})
    set_target_properties(${STATIC_TARGET} PROPERTIES 
        LINK_FLAGS "${EXTRA_LINK_FLAGS}" 
        COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DPSINAD_BUILDING_STATIC_LIBRARY -DPTHREAD_BUILDING_STATIC_LIBRARY")
endif()

if(NOT APPLE AND NOT MSVC)
    find_library(DL_LIBRARY dl)
    if(DL_LIBRARY)
        target_link_libraries(${SHARED_TARGET} PRIVATE ${DL_LIBRARY} ${PTHREADS_LIB})
        target_link_libraries(${STATIC_TARGET} PRIVATE ${DL_LIBRARY} ${PTHREADS_LIB})
        mark_as_advanced(DL_LIBRARY)
    else()
        if (ANDROID AND PSINAD_BUILD_SHARED_LIB)
            target_link_libraries(${SHARED_TARGET} PRIVATE ${PTHREADS_LIB} cpufeatures)
        endif()
    endif()
endif()

if(PSINAD_BUILD_SHARED_LIB)
    install(TARGETS ${SHARED_TARGET} RUNTIME DESTINATION lib)
endif()
if(PSINAD_BUILD_STATIC_LIB)
    install(TARGETS ${STATIC_TARGET} ARCHIVE DESTINATION lib)
endif()
file(GLOB_RECURSE CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
file(GLOB_RECURSE TOP_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/kids/*.h)
install(FILES ${CORE_HEADERS} DESTINATION include)
install(FILES ${TOP_HEADERS} DESTINATION include/kids)

# C++ backends
option(PSINAD_BUILD_CXX_BACKENDS "Build backend for C++" ON)
if(PSINAD_BUILD_CXX_BACKENDS)
    find_package(MPIMod REQUIRED)
    # find_package(MKLMod REQUIRED)

    set(BUILD_TESTING OFF CACHE INTERNAL "Disable unit testing in gflags/glog/catch2")
    set(WITH_GFLAGS OFF CACHE INTERNAL "Option for gflags")
    set(WITH_GTEST OFF CACHE INTERNAL "Option for gtest")

    # On some old computers (e.g., old CentOS), after GCC 5.1, it adds a new realization for
    # std::string & std::list in `libstdc++` by adding `__cxx11` namespace.
    # Some older compilers don't enable this option by default, but it is needed by gflags/glog here.
    add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)

    add_subdirectory(thirdpart/gflags)
    add_subdirectory(thirdpart/glog)
    add_subdirectory(thirdpart/fparser)

    set(BUILD_TESTING ON CACHE INTERNAL "Restore unit testing")
    set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "Disable shared library")

    if(MPIMod_FOUND AND FALSE)
        add_executable(psinad_mpi 
            src/cpprun/psinad_mpi.cpp src/cpprun/mpi_guard.cpp src/cpprun/defined_gflags.cpp)
        target_link_libraries(kids_mpi PRIVATE
            ${SHARED_TARGET}
            gflags::gflags
            glog::glog
            MPIMod::MPI
        )
        install(TARGETS psinad_mpi RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
        
        add_executable(psidyn_mpi 
            src/cpprun/psidyn_mpi.cpp src/cpprun/mpi_guard.cpp src/cpprun/defined_gflags.cpp)
        target_link_libraries(psidyn_mpi PRIVATE
            ${SHARED_TARGET}
            gflags::gflags
            glog::glog
            MPIMod::MPI
        )
        install(TARGETS psidyn_mpi RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
        
        add_executable(psisamp_mpi 
            src/cpprun/psisamp_mpi.cpp src/cpprun/mpi_guard.cpp src/cpprun/defined_gflags.cpp)
        target_link_libraries(psisamp_mpi PRIVATE
            ${SHARED_TARGET}
            gflags::gflags
            glog::glog
            MPIMod::MPI
        )
        install(TARGETS psisamp_mpi RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
    endif()
    add_executable(psinad
        src/cpprun/psinad.cpp src/cpprun/simple_guard.cpp src/cpprun/defined_gflags.cpp)
    target_link_libraries(psinad PRIVATE
        ${SHARED_TARGET}
        gflags::gflags
        glog::glog
    )
    install(TARGETS psinad RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
    
    add_executable(psidyn 
        src/cpprun/psidyn.cpp src/cpprun/simple_guard.cpp src/cpprun/defined_gflags.cpp)
    target_link_libraries(psidyn PRIVATE
        ${SHARED_TARGET}
        gflags::gflags
        glog::glog
    )
    install(TARGETS psidyn RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

    add_executable(psisamp 
        src/cpprun/psisamp.cpp src/cpprun/simple_guard.cpp src/cpprun/defined_gflags.cpp)
    target_link_libraries(psisamp PRIVATE
        ${SHARED_TARGET}
        gflags::gflags
        glog::glog
    )
    install(TARGETS psisamp RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
endif()

# Python wrappers
option(PSINAD_BUILD_PYTHON_WRAPPERS "Build wrappers for Python" ON)
if(PSINAD_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(thirdpart/pybind11)
    add_subdirectory(wrappers/python)
endif()

# Tests
option(PSINAD_BUILD_TESTS "Build tests" OFF)
if(PSINAD_BUILD_TESTS)
    set(CATCH_BUILD_TESTING OFF CACHE INTERNAL "")
    set(CATCH_BUILD_EXAMPLES OFF CACHE INTERNAL "")
    set(CATCH_BUILD_EXTRA_TESTS OFF CACHE INTERNAL "")
    add_subdirectory(thirdpart/catch2)
    enable_testing()
    include_directories(${CMAKE_SOURCE_DIR}/tests)
    add_subdirectory(tests)
endif()

# Examples
option(PSINAD_BUILD_EXAMPLES "Build example executables" OFF)
if(PSINAD_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
