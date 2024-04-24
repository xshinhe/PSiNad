set(KIDS_PYTHON_STAGING_DIR "${CMAKE_BINARY_DIR}/python"
    CACHE PATH "Temporary staging area for Python API wrappers")
mark_as_advanced(KIDS_PYTHON_STAGING_DIR)

# Create package directory structure
file(MAKE_DIRECTORY ${KIDS_PYTHON_STAGING_DIR}/pykids)
file(MAKE_DIRECTORY ${KIDS_PYTHON_STAGING_DIR}/tests)

##############################################################################
### Identify files that need to be copied from source area to staging area ###
##############################################################################
set(STAGING_OUTPUT_FILES "") # Will contain all required package files
file(GLOB STAGING_INPUT_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/MANIFEST.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/README.txt"
)
# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py ${KIDS_PYTHON_STAGING_DIR}/setup.py)

###########################################################
### Check the git revision of the source, and write it  ###
### to a python file in the in the staging directory    ###
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


set(CMAKE_CXX_FLAGS "-Os -fopenmp -march=native -fPIC")

add_library(pykids SHARED libpykids.cpp)

target_link_libraries(
    pykids
    pybind11::embed
    ${SHARED_TARGET}
)


###########################################################################
### Run python setup.py indirectly, so we can set environment variables ###
###########################################################################

# set(PYTHON_SETUP_COMMAND build)
# configure_file(
#     pysetup.cmake.in
#     "${CMAKE_CURRENT_BINARY_DIR}/pysetupbuild.cmake"
#     @ONLY
# )
# add_custom_command(
#     OUTPUT ${KIDS_PYTHON_STAGING_DIR}/build
#     COMMAND ${CMAKE_COMMAND}
#     ARGS -P "${CMAKE_CURRENT_BINARY_DIR}/pysetupbuild.cmake"
#     DEPENDS
#         ${KIDS_PYTHON_STAGING_DIR}/setup.py
#         "${CMAKE_CURRENT_BINARY_DIR}/pysetupbuild.cmake"
#         ${SHARED_TARGET}
#         ${STAGING_OUTPUT_FILES}
#     WORKING_DIRECTORY ${KIDS_PYTHON_STAGING_DIR}
#     COMMENT "Building OpenMM Python binary module..."
# )
# add_custom_target(BuildModule ALL DEPENDS ${STAGING_OUTPUT_FILES})


# set(KIDS_PYTHON_USER_INSTALL OFF CACHE BOOL
#     "Whether to install OpenMM Python binary module into the user site-packages directory")
# mark_as_advanced(KIDS_PYTHON_USER_INSTALL)

# if(KIDS_PYTHON_USER_INSTALL)
#     # Install binary module to user location
#     set(PYTHON_SETUP_COMMAND "install --user")
# else()
#     # Install binary module to system location
#     set(PYTHON_SETUP_COMMAND "install --root=\$ENV{DESTDIR}/")
# endif()
# configure_file(pysetup.cmake.in
#     "${CMAKE_CURRENT_BINARY_DIR}/pysetupinstall.cmake" @ONLY)
# add_custom_target(PythonInstall
#     COMMAND ${CMAKE_COMMAND}
#         -P "${CMAKE_CURRENT_BINARY_DIR}/pysetupinstall.cmake"
#     DEPENDS ${STAGING_OUTPUT_FILES} ${KIDS_PYTHON_STAGING_DIR}/build "${CMAKE_CURRENT_BINARY_DIR}/pysetupinstall.cmake"
#     WORKING_DIRECTORY ${KIDS_PYTHON_STAGING_DIR}
#     COMMENT "Installing OpenMM Python binary module..."
# )
