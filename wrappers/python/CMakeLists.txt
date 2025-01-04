set(PSINAD_PYTHON_STAGING_DIR "${CMAKE_BINARY_DIR}/python"
    CACHE PATH "Temporary path for Python API wrappers")
mark_as_advanced(PSINAD_PYTHON_STAGING_DIR)

# Create package directory structure
file(MAKE_DIRECTORY ${PSINAD_PYTHON_STAGING_DIR}/pyPSiNad)
file(MAKE_DIRECTORY ${PSINAD_PYTHON_STAGING_DIR}/pybind11)
file(MAKE_DIRECTORY ${PSINAD_PYTHON_STAGING_DIR}/swig)
file(MAKE_DIRECTORY ${PSINAD_PYTHON_STAGING_DIR}/ext)
file(MAKE_DIRECTORY ${PSINAD_PYTHON_STAGING_DIR}/plugin)
file(MAKE_DIRECTORY ${PSINAD_PYTHON_STAGING_DIR}/tests)

set(STAGING_OUTPUT_FILES "")
file(GLOB STAGING_INPUT_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/MANIFEST.in"
)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py ${PSINAD_PYTHON_STAGING_DIR}/setup.py)

###########################################################

set(SUBDIRS pyPSiNad pybind11 swig ext plugin tests)
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
    set(outfile "${PSINAD_PYTHON_STAGING_DIR}/${INIT_FILE}")
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
set(SWIG_STAGING_DIR "${PSINAD_PYTHON_STAGING_DIR}/swig")
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
    "${SWIG_STAGING_DIR}/PSINAD_docstring.i"
    "${SWIG_STAGING_DIR}/PSINADSwigHeaders.i"
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
        -d PSINAD_docstring.i
        -o PSINADSwigHeaders.i
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
    "${SWIG_STAGING_DIR}/PSINADSwigHeaders.i" 
    "${SWIG_STAGING_DIR}/lib/typemaps.i")
add_custom_command(
    OUTPUT "${SWIG_STAGING_DIR}/PSINADSwig.cxx" "${PSINAD_PYTHON_STAGING_DIR}/pyPSiNad/pyPSiNad.py"
    COMMAND ${SWIG_EXECUTABLE}
        -python -c++
        -outdir "${PSINAD_PYTHON_STAGING_DIR}/pyPSiNad"
        -o PSINADSwig.cxx
        PSINAD.i
    WORKING_DIRECTORY "${SWIG_STAGING_DIR}"
    DEPENDS
        "${SWIG_STAGING_DIR}/PSINAD.i"
        ${SWIG_INPUT_FILES_REF}
        ${SWIG_INPUT_FILES2}
    COMMENT "Creating Python module sources with swig..."
)
add_custom_target(PSINADSwig DEPENDS
    "${SWIG_STAGING_DIR}/PSINADSwig.cxx"
    "${PSINAD_PYTHON_STAGING_DIR}/pyPSiNad/pyPSiNad.py")
set (STAGING_OUTPUT_FILES ${STAGING_OUTPUT_FILES}
    "${PSINAD_PYTHON_STAGING_DIR}/swig/PSINADSwig.cxx"
    "${PSINAD_PYTHON_STAGING_DIR}/pyPSiNad/pyPSiNad.py")

########################################################
# PYBIND11 STANDALONE
########################################################
# add_subdirectory(pybind11)
# set(CMAKE_CXX_FLAGS "-Os -fopenmp -march=native -fPIC")

add_library(pyPSiNad_v1 SHARED pybind11/libpyPSiNad_v1.cpp)
if(APPLE)
    set(COMPILE_FLAGS_USED "-Os -march=native -fPIC")
else()
    set(COMPILE_FLAGS_USED "-Os -fopenmp -march=native -fPIC")
endif()

set_target_properties(pyPSiNad_v1 PROPERTIES
    LINK_FLAGS "${EXTRA_LINK_FLAGS}"
    COMPILE_FLAGS ${COMPILE_FLAGS_USED}
    LIBRARY_OUTPUT_DIRECTORY "${PSINAD_PYTHON_STAGING_DIR}/pyPSiNad"
    )
target_link_libraries(
    pyPSiNad_v1
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
SET(WRAPPER_INCLUDE_DIRS 
    ${WRAPPER_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/thirdpart
    ${WRAPPER_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/thirdpart/Eigen
)

set(PYTHON_SETUP_COMMAND build)
configure_file(pysetup.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/pysetupbuild.cmake" @ONLY)
set(PYTHON_SETUP_COMMAND bdist)
configure_file(pysetup.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/pysetupbdist.cmake" @ONLY)
set(PYTHON_SETUP_COMMAND sdist)
configure_file(pysetup.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/pysetupsdist.cmake" @ONLY)
set(PSINAD_PYTHON_USER_INSTALL OFF CACHE BOOL
    "Whether to install Python binary module into the user site-packages directory")
mark_as_advanced(PSINAD_PYTHON_USER_INSTALL)
if(PSINAD_PYTHON_USER_INSTALL)
    set(PYTHON_SETUP_COMMAND "install --user")
else()
    set(PYTHON_SETUP_COMMAND "install --root=\$ENV{DESTDIR}/")
endif()
configure_file(pysetup.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/pysetupinstall.cmake" @ONLY)
configure_file(pipsetup.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/pipsetupinstall.cmake" @ONLY)

add_custom_command(
    OUTPUT ${PSINAD_PYTHON_STAGING_DIR}/build
    COMMAND ${CMAKE_COMMAND}
    ARGS -P "${CMAKE_CURRENT_BINARY_DIR}/pysetupbuild.cmake"
    DEPENDS
        ${PSINAD_PYTHON_STAGING_DIR}/setup.py
        "${CMAKE_CURRENT_BINARY_DIR}/pysetupbuild.cmake"
        ${SHARED_TARGET}
        ${STAGING_OUTPUT_FILES}
    WORKING_DIRECTORY ${PSINAD_PYTHON_STAGING_DIR}
    COMMENT "Building PSINAD Python binary module..."
)
add_custom_target(BuildModule ALL DEPENDS ${STAGING_OUTPUT_FILES})
add_custom_target(PythonBdist
    COMMAND ${CMAKE_COMMAND}
        -P "${CMAKE_CURRENT_BINARY_DIR}/pysetupbdist.cmake"
    DEPENDS ${STAGING_OUTPUT_FILES} "${CMAKE_CURRENT_BINARY_DIR}/pysetupbdist.cmake"
    WORKING_DIRECTORY ${PSINAD_PYTHON_STAGING_DIR}
    COMMENT "Packaging binary distribution package (bdist)..."
)
add_custom_target(PythonSdist
    COMMAND ${CMAKE_COMMAND}
        -P "${CMAKE_CURRENT_BINARY_DIR}/pysetupsdist.cmake"
    DEPENDS ${STAGING_OUTPUT_FILES} "${CMAKE_CURRENT_BINARY_DIR}/pysetupsdist.cmake"
    WORKING_DIRECTORY ${PSINAD_PYTHON_STAGING_DIR}
    COMMENT "Packaging source distribution package (sdist)..."
)
add_custom_target(PythonInstall
    COMMAND ${CMAKE_COMMAND}
        -P "${CMAKE_CURRENT_BINARY_DIR}/pysetupinstall.cmake"
    DEPENDS ${STAGING_OUTPUT_FILES} ${PSINAD_PYTHON_STAGING_DIR}/build "${CMAKE_CURRENT_BINARY_DIR}/pysetupinstall.cmake"
    WORKING_DIRECTORY ${PSINAD_PYTHON_STAGING_DIR}
    COMMENT "Installing PSINAD Python binary module..."
)