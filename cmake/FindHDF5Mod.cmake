#[==[

FindHDF5Mod
-----------

Dependencies
^^^^^^^^^^^^
`HDF5_ROOT`: the variable containing the path of the library root directory,
which must have
    bin/<hdf5 CLIs>
    include/hdf5.h
    lib/<hdf5 libraries> (e.g., libhdf5.so)

Resulted variables
^^^^^^^^^^^^^^^^^^
HDF5Mod_FOUND

Cache entry
^^^^^^^^^^^
HDF5Mod_INCLUDE_DIR
HDF5Mod_LIBRARY

Imported library
^^^^^^^^^^^^^^^^
HDF5Mod::HDF5

#]==]

set(hdf5_mod_find_flags)
if(HDF5Mod_FIND_REQUIRED)
    list(APPEND hdf5_mod_find_flags "REQUIRED")
endif()
if(HDF5Mod_FIND_QUIETLY)
    list(APPEND hdf5_mod_find_flags "QUIET")
endif()

find_package(HDF5 COMPONENTS C HL NO_MODULE) # MODIFIED TO FIX PROBLEM OF offline installed HDF5

set(HDF5Mod_FOUND ${HDF5_FOUND} CACHE BOOL "Whether or not HDF5 Module is found?")

if(HDF5Mod_FOUND AND NOT TARGET HDF5Mod::HDF5)
    message(STATUS "Imported library HDF5Mod::HDF5")
    message("   Include dirs: ${HDF5_INCLUDE_DIRS}")
    message("   Libraries: ${HDF5_LIBRARIES}")

    set(HDF5Mod_INCLUDE_DIR "${HDF5_INCLUDE_DIRS}" CACHE STRING "HDF5 include directories")
    set(HDF5Mod_LIBRARY "${HDF5_LIBRARIES}" CACHE STRING "HDF5 libraries")

    add_library(HDF5Mod::HDF5 INTERFACE IMPORTED)
    target_include_directories(HDF5Mod::HDF5
        INTERFACE ${HDF5_INCLUDE_DIRS}
    )
    target_link_libraries(HDF5Mod::HDF5
        INTERFACE ${HDF5_LIBRARIES}
    )

elseif(DEFINED CUSTOMIZED_HDF5_DIR AND NOT TARGET HDF5Mod::HDF5)
    find_path(HDF5Mod_INCLUDE_PATH NAMES "hdf5.h"
        PATHS ${CUSTOMIZED_HDF5_DIR}
        PATH_SUFFIXES
        include
    )
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(HDF5Mod
        DEFAULT_MSG
        HDF5Mod_INCLUDE_PATH
    )

    set(HDF5Mod_LIBRARIES
            libhdf5.so       # NOTE: it should be parallel version
    )

    add_library(HDF5Mod::HDF5 INTERFACE IMPORTED)
    target_include_directories(HDF5Mod::HDF5
        INTERFACE "${HDF5Mod_INCLUDE_PATH}"
    )
    target_link_directories(HDF5Mod::HDF5
        INTERFACE "${CUSTOMIZED_HDF5_DIR}"
        INTERFACE "${CUSTOMIZED_HDF5_DIR}/lib"
    )
    target_link_libraries(HDF5Mod::HDF5
        INTERFACE "${HDF5Mod_LIBRARIES}"
    )
endif()
