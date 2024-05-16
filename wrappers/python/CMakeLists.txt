set(KIDS_PYTHON_STAGING_DIR "${CMAKE_BINARY_DIR}/python"
    CACHE PATH "Temporary path for Python API wrappers")
mark_as_advanced(KIDS_PYTHON_STAGING_DIR)

# Create package directory structure
file(MAKE_DIRECTORY ${KIDS_PYTHON_STAGING_DIR}/pykids)
file(MAKE_DIRECTORY ${KIDS_PYTHON_STAGING_DIR}/pybind11)
file(MAKE_DIRECTORY ${KIDS_PYTHON_STAGING_DIR}/swig)
file(MAKE_DIRECTORY ${KIDS_PYTHON_STAGING_DIR}/ext)
file(MAKE_DIRECTORY ${KIDS_PYTHON_STAGING_DIR}/plugin)
file(MAKE_DIRECTORY ${KIDS_PYTHON_STAGING_DIR}/tests)

set(STAGING_OUTPUT_FILES "")
file(GLOB STAGING_INPUT_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/MANIFEST.in"
)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py ${KIDS_PYTHON_STAGING_DIR}/setup.py)

###########################################################
execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE rev_hash_str
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
if(NOT rev_hash_str)
  set(rev_hash_str "Unknown")
endif()
file(WRITE "${KIDS_PYTHON_STAGING_DIR}/pykids/version.py" "git_revision = '${rev_hash_str}'\n")

set(SUBDIRS pykids pybind11 swig ext plugin tests)
foreach(SUBDIR ${SUBDIRS})
    file(GLOB_RECURSE STAGING_INPUT_FILES1 RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*README.txt"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.pyx"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.pxd"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.cc"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.i"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.ini"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.sh"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.xml"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.json"
    )
    foreach(file ${STAGING_INPUT_FILES1})
        set(STAGING_INPUT_FILES ${STAGING_INPUT_FILES} "${file}")
    endforeach(file ${STAGING_INPUT_FILES1})
endforeach(SUBDIR ${SUBDIRS})
foreach(INIT_FILE ${STAGING_INPUT_FILES})
    set(infile "${CMAKE_CURRENT_SOURCE_DIR}/${INIT_FILE}")
    set(outfile "${KIDS_PYTHON_STAGING_DIR}/${INIT_FILE}")
    add_custom_command(
        OUTPUT "${outfile}"
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${infile}" "${outfile}"
        DEPENDS "${infile}"
        COMMENT "CMake-copying file ${infile} to ${outfile}")
    set(STAGING_OUTPUT_FILES ${STAGING_OUTPUT_FILES} "${outfile}")
endforeach(INIT_FILE ${STAGING_INPUT_FILES})


########################################################
# SWIG/DOXYGEN
########################################################
find_program(SWIG_EXECUTABLE swig PATHS
    "$ENV{HOME}/bin" "/Users/builder/bin" "/home/builder/bin"
    NO_DEFAULT_PATH)
find_program(SWIG_EXECUTABLE swig)
MARK_AS_ADVANCED(SWIG_VERSION)
if(SWIG_EXECUTABLE)
    execute_process(COMMAND ${SWIG_EXECUTABLE} -version
        OUTPUT_VARIABLE SWIG_VERSION_output ERROR_VARIABLE SWIG_VERSION_output)
    string(REGEX REPLACE ".*SWIG Version[^0-9.]*\([0-9.]+\).*" "\\1"
        SWIG_VERSION_output "${SWIG_VERSION_output}")
    set(SWIG_VERSION ${SWIG_VERSION_output} CACHE STRING "Swig version" FORCE)
else(SWIG_EXECUTABLE)
    set(SWIG_VERSION "0.0.0" CACHE STRING "Swig version" FORCE)
endif(SWIG_EXECUTABLE)
if(SWIG_VERSION VERSION_LESS "3.0.5")
    message(SEND_ERROR "Swig version must be 3.0.5 or greater! (You have ${SWIG_VERSION})")
endif(SWIG_VERSION VERSION_LESS "3.0.5")
find_package(Doxygen REQUIRED)
mark_as_advanced(CLEAR DOXYGEN_EXECUTABLE)
set(SWIG_STAGING_DIR "${KIDS_PYTHON_STAGING_DIR}/swig")
file(MAKE_DIRECTORY ${SWIG_STAGING_DIR}/lib)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/swig/doxygen/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/swig/doxygen/Doxyfile
    @ONLY
)
add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/swig/doxygen/xml/index.xml"
    COMMAND "${DOXYGEN_EXECUTABLE}"
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/swig/doxygen/Doxyfile"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/swig/doxygen"
    COMMENT "Parsing header files with doxygen..."
)
set(SWIG_INPUT_FILES_REF
    "${SWIG_STAGING_DIR}/lib/pythonprepend.i"
    "${SWIG_STAGING_DIR}/lib/pythonappend.i"
    "${SWIG_STAGING_DIR}/KIDS_docstring.i"
    "${SWIG_STAGING_DIR}/KIDSSwigHeaders.i"
)
set(SWIG_INPUT_FILES2
    "${SWIG_STAGING_DIR}/lib/exceptions.i"
    "${SWIG_STAGING_DIR}/lib/extend.i"
    "${SWIG_STAGING_DIR}/lib/features.i"
    "${SWIG_STAGING_DIR}/lib/header.i"
    "${SWIG_STAGING_DIR}/lib/pythoncode.i"
    "${SWIG_STAGING_DIR}/lib/typemaps.i"
)
add_custom_command(
    OUTPUT ${SWIG_INPUT_FILES_REF}
    COMMAND ${PYTHON_EXECUTABLE} "${SWIG_STAGING_DIR}/swigInputBuilder.py"
        -i "${CMAKE_CURRENT_BINARY_DIR}/swig/doxygen/xml"
        -c swigInputConfig.py
        -d KIDS_docstring.i
        -o KIDSSwigHeaders.i
        -a lib/pythonprepend.i
        -z lib/pythonappend.i
        -v "${SWIG_VERSION}"
    WORKING_DIRECTORY "${SWIG_STAGING_DIR}"
    DEPENDS
        "${SWIG_STAGING_DIR}/swigInputConfig.py"
        "${SWIG_STAGING_DIR}/swigInputBuilder.py"
        "${CMAKE_CURRENT_BINARY_DIR}/swig/doxygen/xml/index.xml"
    COMMENT "Creating Python swig input files..."
)
install_files(/include/swig FILES 
    "${SWIG_STAGING_DIR}/KIDSSwigHeaders.i" 
    "${SWIG_STAGING_DIR}/lib/typemaps.i")
add_custom_command(
    OUTPUT "${SWIG_STAGING_DIR}/KIDSSwig.cxx" "${KIDS_PYTHON_STAGING_DIR}/pykids/pykids.py"
    COMMAND ${SWIG_EXECUTABLE}
        -python -c++
        -outdir "${KIDS_PYTHON_STAGING_DIR}/pykids"
        -o KIDSSwig.cxx
        KIDS.i
    WORKING_DIRECTORY "${SWIG_STAGING_DIR}"
    DEPENDS
        "${SWIG_STAGING_DIR}/KIDS.i"
        ${SWIG_INPUT_FILES_REF}
        ${SWIG_INPUT_FILES2}
    COMMENT "Creating Python module sources with swig..."
)
add_custom_target(KIDSSwig DEPENDS
    "${SWIG_STAGING_DIR}/KIDSSwig.cxx"
    "${KIDS_PYTHON_STAGING_DIR}/pykids/pykids.py")
set (STAGING_OUTPUT_FILES ${STAGING_OUTPUT_FILES}
    "${KIDS_PYTHON_STAGING_DIR}/swig/KIDSSwig.cxx"
    "${KIDS_PYTHON_STAGING_DIR}/pykids/pykids.py")

########################################################
# PYBIND11 STANDALONE
########################################################
# add_subdirectory(pybind11)
# set(CMAKE_CXX_FLAGS "-Os -fopenmp -march=native -fPIC")

add_library(pykids_v1 SHARED pybind11/libpykids_v1.cpp)
set_target_properties(pykids_v1 PROPERTIES
    LINK_FLAGS "${EXTRA_LINK_FLAGS}"
    COMPILE_FLAGS "-Os -fopenmp -march=native -fPIC"
    LIBRARY_OUTPUT_DIRECTORY "${KIDS_PYTHON_STAGING_DIR}/pykids"
    )
target_link_libraries(
    pykids_v1
    pybind11::embed
    ${SHARED_TARGET}
)

###########################################################################
### Run python setup.py indirectly, so we can set environment variables ###
###########################################################################

SET(WRAPPER_BASE_SUBDIRS src/core src/kernels src/models src/solvers)
SET(WRAPPER_INCLUDE_DIRS)
FOREACH(subdir ${WRAPPER_BASE_SUBDIRS})
    SET(WRAPPER_INCLUDE_DIRS ${WRAPPER_INCLUDE_DIRS}
                         ${CMAKE_SOURCE_DIR}/${subdir}/include
                         ${CMAKE_SOURCE_DIR}/${subdir}/include/kids)
ENDFOREACH(subdir)
SET(WRAPPER_INCLUDE_DIRS ${WRAPPER_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/thirdpart)

set(PYTHON_SETUP_COMMAND build)
configure_file(pysetup.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/pysetupbuild.cmake" @ONLY)
set(PYTHON_SETUP_COMMAND bdist)
configure_file(pysetup.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/pysetupbdist.cmake" @ONLY)
set(PYTHON_SETUP_COMMAND sdist)
configure_file(pysetup.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/pysetupsdist.cmake" @ONLY)
set(KIDS_PYTHON_USER_INSTALL OFF CACHE BOOL
    "Whether to install Python binary module into the user site-packages directory")
mark_as_advanced(KIDS_PYTHON_USER_INSTALL)
if(KIDS_PYTHON_USER_INSTALL)
    set(PYTHON_SETUP_COMMAND "install --user")
else()
    set(PYTHON_SETUP_COMMAND "install --root=\$ENV{DESTDIR}/")
endif()
configure_file(pysetup.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/pysetupinstall.cmake" @ONLY)

add_custom_command(
    OUTPUT ${KIDS_PYTHON_STAGING_DIR}/build
    COMMAND ${CMAKE_COMMAND}
    ARGS -P "${CMAKE_CURRENT_BINARY_DIR}/pysetupbuild.cmake"
    DEPENDS
        ${KIDS_PYTHON_STAGING_DIR}/setup.py
        "${CMAKE_CURRENT_BINARY_DIR}/pysetupbuild.cmake"
        ${SHARED_TARGET}
        ${STAGING_OUTPUT_FILES}
    WORKING_DIRECTORY ${KIDS_PYTHON_STAGING_DIR}
    COMMENT "Building KIDS Python binary module..."
)
add_custom_target(BuildModule ALL DEPENDS ${STAGING_OUTPUT_FILES})
add_custom_target(PythonBdist
    COMMAND ${CMAKE_COMMAND}
        -P "${CMAKE_CURRENT_BINARY_DIR}/pysetupbdist.cmake"
    DEPENDS ${STAGING_OUTPUT_FILES} "${CMAKE_CURRENT_BINARY_DIR}/pysetupbdist.cmake"
    WORKING_DIRECTORY ${KIDS_PYTHON_STAGING_DIR}
    COMMENT "Packaging binary distribution package (bdist)..."
)
add_custom_target(PythonSdist
    COMMAND ${CMAKE_COMMAND}
        -P "${CMAKE_CURRENT_BINARY_DIR}/pysetupsdist.cmake"
    DEPENDS ${STAGING_OUTPUT_FILES} "${CMAKE_CURRENT_BINARY_DIR}/pysetupsdist.cmake"
    WORKING_DIRECTORY ${KIDS_PYTHON_STAGING_DIR}
    COMMENT "Packaging source distribution package (sdist)..."
)
add_custom_target(PythonInstall
    COMMAND ${CMAKE_COMMAND}
        -P "${CMAKE_CURRENT_BINARY_DIR}/pysetupinstall.cmake"
    DEPENDS ${STAGING_OUTPUT_FILES} ${KIDS_PYTHON_STAGING_DIR}/build "${CMAKE_CURRENT_BINARY_DIR}/pysetupinstall.cmake"
    WORKING_DIRECTORY ${KIDS_PYTHON_STAGING_DIR}
    COMMENT "Installing KIDS Python binary module..."
)