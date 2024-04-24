cmake_minimum_required(VERSION 3.17)

cmake_policy(SET CMP0042 OLD)
cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0005 NEW)
cmake_policy(SET CMP0011 NEW)

# You should choose a different directory from source for compilation
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "
        FATAL: In-source builds are not allowed. 
        You should create a separate directory for build files.
    ")
endif()
set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL  "Whether to verbose in compilation")

# Don't create a new project name if this is part of the parent directory
if(NOT PROJECT_NAME)
    project(kids LANGUAGES C CXX Fortran)
endif(NOT PROJECT_NAME)
enable_language(Fortran)

find_package(PythonInterp)
mark_as_advanced(CLEAR PYTHON_EXECUTABLE)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/thirdpart/catch2/contrib)

# Determine what architecture we are compiling for.
include(TargetArch)
target_architecture(TARGET_ARCH)
if ("${TARGET_ARCH}" MATCHES "x86_64|i386")
    set(X86 ON)
endif()
if ("${TARGET_ARCH}" MATCHES "arm")
    set(ARM ON)
    if ("${TARGET_ARCH}" MATCHES "armv8")
        add_definitions(-D__ARM64__=1)
    else()
        add_definitions(-D__ARM__=1)
    endif()
endif()
if ("${TARGET_ARCH}" MATCHES "ppc")
    set(PPC ON)
    add_definitions(-D__PPC__=1)
endif()


# What type of libraries to build
set(KIDS_BUILD_SHARED_LIB ON CACHE BOOL  "Whether to build shared KIDS libraries")
set(KIDS_BUILD_STATIC_LIB ON CACHE BOOL  "Whether to build static KIDS libraries")

# processing modular configuration (solver / model)
# include(module)
# module_config()

# help to optimize and measure the compilation time
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    message(STATUS "set up ccache ...")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif()

find_package(compile-time-perf) # https://github.com/jrmadsen/compile-time-perf
if(compile-time-perf_FOUND)
    enable_compile_time_perf(kids-ctp)
endif()

include_directories(${CMAKE_SOURCE_DIR}/generate)
include_directories(${CMAKE_SOURCE_DIR}/thirdpart)

################################################################################
# set esssential variables
################################################################################

# set(CMAKE_BUILD_TYPE <Config>)
## Please -DCMAKE_BUILD_TYPE=<Config> : Debug, Release, RelWithDebInfo, MinsizeRel
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -march=native -fPIC")

# ################################################################################
# # configure libraries
# ################################################################################
# set(CUSTOMIZED_MKL_DIR "/opt/intel/oneapi/mkl/latest/")
# set(CUSTOMIZED_MPI_DIR "/opt/intel/oneapi/mpi/latest/")
# set(CUSTOMIZED_HDF5_DIR "/usr/local/HDF_Group/HDF5/1.13.2")

#find_package(MKLMod)

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

# add_subdirectory(core)
# add_subdirectory(kernels)
# add_subdirectory(models)

# add_subdirectory(solvers)
# file(GLOB_RECURSE kids_solver_srcs solvers/*.cpp)
# add_library(cppkids SHARED ${kids_solver_srcs})
# target_link_libraries(cppkids
    # fparser
    # core
    # ker
    # model
# )

####################################

# Where to install KIDS
if(NOT KIDS_INSTALL_PREFIX)
if(WIN32)
    set(KIDS_INSTALL_PREFIX "$ENV{ProgramFiles}/KIDS")
else(WIN32)
    set(KIDS_INSTALL_PREFIX "/usr/local/kids")
endif(WIN32)
endif(NOT KIDS_INSTALL_PREFIX)
if(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})
  set(CMAKE_INSTALL_PREFIX ${KIDS_INSTALL_PREFIX} CACHE PATH "Where to install KIDS" FORCE)
endif(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})

# The source is organized into subdirectories, but we handle them all togethor
# rather than visit them as SUBDIRS.
set(KIDS_SOURCE_SUBDIRS
    src/core
    src/kernels
    src/models
    src/solvers
)
file(GLOB third_src_files  ${CMAKE_CURRENT_SOURCE_DIR}/thirdpart/fparser/*.cpp)

# if(X86 OR ARM)
#     set(KIDS_SOURCE_SUBDIRS ${KIDS_SOURCE_SUBDIRS} libraries/vecmath)
# endif()
# if(WIN32)
#     add_definitions(-DWIN32)
#     set(KIDS_SOURCE_SUBDIRS ${KIDS_SOURCE_SUBDIRS} libraries/pthreads)
# else(WIN32)
#     if (NOT ANDROID)
#         set(PTHREADS_LIB pthread)
#         set(PTHREADS_LIB_STATIC pthread)
#     endif (NOT ANDROID)
# endif(WIN32)

if(ANDROID) # Include CPU-Features for Android
    include_directories(${ANDROID_NDK}/sources/cpufeatures)
endif(ANDROID)

if(APPLE)
    if (NOT CMAKE_OSX_DEPLOYMENT_TARGET)
        set (CMAKE_OSX_DEPLOYMENT_TARGET "10.7" CACHE STRING "The minimum version of OS X to support" FORCE)
    endif (NOT CMAKE_OSX_DEPLOYMENT_TARGET)
    if (NOT CMAKE_OSX_ARCHITECTURES)
        set (CMAKE_OSX_ARCHITECTURES "${CMAKE_SYSTEM_PROCESSOR}" CACHE STRING "The processor architectures to build for" FORCE)
    endif (NOT CMAKE_OSX_ARCHITECTURES)
    if (NOT CMAKE_OSX_SYSROOT)
        execute_process(COMMAND "xcrun" "--show-sdk-path" OUTPUT_VARIABLE XCRUN_OSX_SYSROOT RESULT_VARIABLE XCRUN_OSX_SYSROOT_STATUS OUTPUT_STRIP_TRAILING_WHITESPACE)
        if (XCRUN_OSX_SYSROOT_STATUS EQUAL 0)
            set (CMAKE_OSX_SYSROOT "${XCRUN_OSX_SYSROOT}" CACHE STRING "SDK Path" FORCE)
        endif (XCRUN_OSX_SYSROOT_STATUS EQUAL 0)
        unset(XCRUN_OSX_SYSROOT)
        unset(XCRUN_OSX_SYSROOT_STATUS)
    endif (NOT CMAKE_OSX_SYSROOT)

    # Improve the linking behavior of Mac libraries
    set (CMAKE_INSTALL_NAME_DIR "@rpath")
    if (X86)
      set(EXTRA_COMPILE_FLAGS "-msse2")
    endif()
    set(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -stdlib=libc++")
elseif(MSVC)
    set(EXTRA_COMPILE_FLAGS)
    # Use warning level 2, not whatever warning level CMake picked.
    string(REGEX REPLACE "/W[0-4]" "/W2" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    # Explicitly suppress warnings 4305 and 4244.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4305 /wd4244")
elseif(X86)
    set(EXTRA_COMPILE_FLAGS "-msse2")
else()
    set(EXTRA_COMPILE_FLAGS)
endif()

if(UNIX AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release build" FORCE)
endif (UNIX AND NOT CMAKE_BUILD_TYPE)

if (NOT CMAKE_CXX_FLAGS_DEBUG)
    set(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "To use when CMAKE_BUILD_TYPE=Debug" FORCE)
endif (NOT CMAKE_CXX_FLAGS_DEBUG)

if (NOT CMAKE_CXX_FLAGS_RELEASE)
    if (ANDROID)
        set(CMAKE_CXX_FLAGS_RELEASE "-mfloat-abi=softfp -march=armv7-a -mfpu=neon -funsafe-math-optimizations -O3 -DNDEBUG" CACHE STRING
        "To use when CMAKE_BUILD_TYPE=Release" FORCE)
    else (ANDROID)
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING
        "To use when CMAKE_BUILD_TYPE=Release" FORCE)
    endif (ANDROID)
endif (NOT CMAKE_CXX_FLAGS_RELEASE)



# Collect up information about the version of the KIDS library we're building
# and make it available to the code so it can be built into the binaries.

# Processing repositary fingerprint
include(repo)
get_repo_info(GIT_BRANCH GIT_DATE GIT_HASH GIT_STATUS)
configure_file(
  ${CMAKE_SOURCE_DIR}/version/version.h.in
  ${CMAKE_SOURCE_DIR}/generate/version.h
  @ONLY
)

set(KIDS_LIBRARY_NAME KIDS)
set(KIDS_MAJOR_VERSION 0)
set(KIDS_MINOR_VERSION 1)
set(KIDS_BUILD_VERSION 1)

add_definitions(-DKIDS_LIBRARY_NAME=${KIDS_LIBRARY_NAME}
                -DKIDS_MAJOR_VERSION=${KIDS_MAJOR_VERSION}
                -DKIDS_MINOR_VERSION=${KIDS_MINOR_VERSION}
                -DKIDS_BUILD_VERSION=${KIDS_BUILD_VERSION})

# -DKIDS_TYPE has to be defined in the target subdirectories.
# -Dsimbody_EXPORTS defined automatically when Windows DLL build is being done.

# Report the version number to the CMake UI
set(KIDS_VERSION
    "${KIDS_MAJOR_VERSION}.${KIDS_MINOR_VERSION}.${KIDS_BUILD_VERSION}"
     CACHE STRING "This is the version of KIDS which will be built." FORCE)
mark_as_advanced(KIDS_VERSION)

set(SHARED_TARGET ${KIDS_LIBRARY_NAME})
set(STATIC_TARGET ${KIDS_LIBRARY_NAME}_static)

## If no one says otherwise, change the executable path to drop into the same binary
## location as the DLLs so that the test cases will use the just-build DLLs.
if(NOT EXECUTABLE_OUTPUT_PATH)
  set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}
      CACHE INTERNAL "Single output directory for building all executables.")
endif(NOT EXECUTABLE_OUTPUT_PATH)
if(NOT LIBRARY_OUTPUT_PATH)
  set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}
      CACHE INTERNAL "Single output directory for building all libraries.")
endif(NOT LIBRARY_OUTPUT_PATH)
set(${PROJECT_NAME}_EXECUTABLE_DIR ${EXECUTABLE_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}) #??
set(${PROJECT_NAME}_LIBRARY_DIR    ${LIBRARY_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}) #??


# These are all the places to search for header files which are
# to be part of the API.
set(API_INCLUDE_DIRS) # start empty
foreach(subdir ${KIDS_SOURCE_SUBDIRS})
    # append
    set(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}
                         ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include
                         ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include/kids)
endforeach(subdir)

# We'll need both *relative* path names, starting with their API_INCLUDE_DIRS,
# and absolute pathnames.
set(API_REL_INCLUDE_FILES)   # start these out empty
set(API_ABS_INCLUDE_FILES)

foreach(dir ${API_INCLUDE_DIRS})
    file(GLOB fullpaths ${dir}/*.h) # returns full pathnames
    set(API_ABS_INCLUDE_FILES ${API_ABS_INCLUDE_FILES} ${fullpaths})

    foreach(pathname ${fullpaths})
        get_filename_component(filename ${pathname} NAME)
        set(API_REL_INCLUDE_FILES ${API_REL_INCLUDE_FILES} ${dir}/${filename})
    endforeach(pathname)
endforeach(dir)

# collect up source files
set(SOURCE_FILES) # empty
set(SOURCE_INCLUDE_FILES)

foreach(subdir ${KIDS_SOURCE_SUBDIRS})
    file(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*/*.cpp)
    file(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    set(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    set(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in KIDS/include if
    ## KIDS was previously installed there.
    include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
endforeach(subdir)



if(X86)
    set_source_files_properties(${CMAKE_SOURCE_DIR}/libraries/sfmt/src/SFMT.cpp PROPERTIES COMPILE_FLAGS "-DHAVE_SSE2=1")
else()
    set_source_files_properties(${CMAKE_SOURCE_DIR}/libraries/sfmt/src/SFMT.cpp PROPERTIES COMPILE_FLAGS "-UHAVE_SSE2")
endif()
if((X86 OR ARM) AND NOT (WIN32 AND KIDS_BUILD_STATIC_LIB))
    file(GLOB src_files ${CMAKE_CURRENT_SOURCE_DIR}/libraries/asmjit/asmjit/*/*.cpp)
    file(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/libraries/asmjit/*.h)
    set(SOURCE_FILES ${SOURCE_FILES} ${src_files})
    set(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
    include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/libraries/asmjit")
    set(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DLEPTON_USE_JIT")
endif()


include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(EXTRA_LINK_FLAGS ${EXTRA_COMPILE_FLAGS})
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    if (NOT ANDROID)
        set(EXTRA_LINK_FLAGS "${EXTRA_LINK_FLAGS} -Wl,--no-as-needed -lrt")
    endif (NOT ANDROID)
endif (CMAKE_SYSTEM_NAME MATCHES "Linux")
if (MSVC)
    set(EXTRA_LINK_FLAGS)
endif (MSVC)

set(SOURCE_FILES ${SOURCE_FILES} ${third_src_files})

if(KIDS_BUILD_SHARED_LIB)
    add_library(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})
    set_target_properties(${SHARED_TARGET} PROPERTIES LINK_FLAGS "${EXTRA_LINK_FLAGS}" COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DKIDS_BUILDING_SHARED_LIBRARY -DLEPTON_BUILDING_SHARED_LIBRARY -DPTHREAD_BUILDING_SHARED_LIBRARY" SOVERSION "${KIDS_MAJOR_VERSION}.${KIDS_MINOR_VERSION}")
endif(KIDS_BUILD_SHARED_LIB)

if(KIDS_BUILD_STATIC_LIB)
    add_library(${STATIC_TARGET} STATIC ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})
    set(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DKIDS_USE_STATIC_LIBRARIES -DLEPTON_USE_STATIC_LIBRARIES -DPTW32_STATIC_LIB")
    set_target_properties(${STATIC_TARGET} PROPERTIES LINK_FLAGS "${EXTRA_LINK_FLAGS}" COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DKIDS_BUILDING_STATIC_LIBRARY -DLEPTON_BUILDING_STATIC_LIBRARY -DPTHREAD_BUILDING_STATIC_LIBRARY")
endif(KIDS_BUILD_STATIC_LIB)

# if(KIDS_BUILD_C_AND_FORTRAN_WRAPPERS)
#     add_dependencies(${SHARED_TARGET} ApiWrappers)
#     if(KIDS_BUILD_STATIC_LIB)
#         add_dependencies(${STATIC_TARGET} ApiWrappers)
#     endif(KIDS_BUILD_STATIC_LIB)
# endif(KIDS_BUILD_C_AND_FORTRAN_WRAPPERS)

# On Linux need to link to libdl
if(NOT APPLE)
    FIND_LIBRARY(DL_LIBRARY dl)
    if(DL_LIBRARY)
        if(KIDS_BUILD_SHARED_LIB)
            target_link_libraries(${SHARED_TARGET} ${DL_LIBRARY} ${PTHREADS_LIB})
        endif(KIDS_BUILD_SHARED_LIB)
        if(KIDS_BUILD_STATIC_LIB)
            target_link_libraries(${STATIC_TARGET} ${DL_LIBRARY} ${PTHREADS_LIB})
        endif(KIDS_BUILD_STATIC_LIB)
        mark_as_advanced(DL_LIBRARY)
    else(DL_LIBRARY)
        if(KIDS_BUILD_SHARED_LIB)
            if (ANDROID)
                target_link_libraries(${SHARED_TARGET} ${PTHREADS_LIB} cpufeatures)
            else (ANDROID)
                target_link_libraries(${SHARED_TARGET} ${PTHREADS_LIB})
            endif (ANDROID)
        endif(KIDS_BUILD_SHARED_LIB)
    endif(DL_LIBRARY)
endif()

if(KIDS_BUILD_SHARED_LIB)
    install_targets(/lib RUNTIME_DIRECTORY /lib ${SHARED_TARGET})
endif(KIDS_BUILD_SHARED_LIB)
if(KIDS_BUILD_STATIC_LIB)
    install_targets(/lib RUNTIME_DIRECTORY /lib ${STATIC_TARGET})
endif(KIDS_BUILD_STATIC_LIB)
file(GLOB CORE_HEADERS      include/*.h                 */include/*.h)
file(GLOB TOP_HEADERS       include/kids/*.h            */include/kids/*.h)
install_files(/include                  FILES ${CORE_HEADERS})
install_files(/include/kids             FILES ${TOP_HEADERS})

# If API wrappers are being generated, and add them to the build.
# set(KIDS_BUILD_C_AND_FORTRAN_WRAPPERS ON CACHE BOOL "Build wrappers for C and Fortran")
# if(KIDS_BUILD_C_AND_FORTRAN_WRAPPERS)
#     add_subdirectory(wrappers)
#     set(SOURCE_FILES ${SOURCE_FILES} wrappers/KIDSCWrapper.cpp  wrappers/KIDSFortranWrapper.cpp)
#     set_source_files_properties(wrappers/KIDSCWrapper.cpp wrappers/KIDSFortranWrapper.cpp PROPERTIES GENERATED TRUE)
# endif(KIDS_BUILD_C_AND_FORTRAN_WRAPPERS)

# CXX backends
set(KIDS_BUILD_CXX_BACKENDS ON CACHE BOOL "Build backend for C++")
if(KIDS_BUILD_CXX_BACKENDS)
    find_package(MPIMod REQUIRED)

    # add_subdirectory(thirdpart/gflags)
    # add_subdirectory(thirdpart/glog)

    add_executable(kids src/kidsrun/main.cc src/kidsrun/Handler.cpp src/kidsrun/mpi_utils.cpp)

    target_link_libraries(kids
        ${SHARED_TARGET}
        gflags::gflags
        glog::glog
        # MKLMod::MKL
        MPIMod::MPI
    )
endif(KIDS_BUILD_CXX_BACKENDS)

# Python wrappers
set(KIDS_BUILD_PYTHON_WRAPPERS ON CACHE BOOL "Build wrappers for Python")
if(KIDS_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(wrappers/pybind11)
endif(KIDS_BUILD_PYTHON_WRAPPERS)


if (EXECUTABLE_OUTPUT_PATH)
  set (TEST_PATH ${EXECUTABLE_OUTPUT_PATH})
else (EXECUTABLE_OUTPUT_PATH)
  set (TEST_PATH .)
endif (EXECUTABLE_OUTPUT_PATH)

# add_subdirectory(docs)
set(KIDS_BUILD_TESTS OFF CACHE BOOL "Build tests")
if(KIDS_BUILD_TESTS)
    add_subdirectory(thirdpart/catch2)

    enable_testing()
    include_directories(${CMAKE_SOURCE_DIR}/tests)
    add_subdirectory(tests)
endif(KIDS_BUILD_TESTS)

set(KIDS_BUILD_EXAMPLES OFF CACHE BOOL "Build example executables")
if(KIDS_BUILD_EXAMPLES)
  # add_subdirectory(examples)
endif(KIDS_BUILD_EXAMPLES)
