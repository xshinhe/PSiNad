#[==[

FindMPIMod
-----------

Resulted variables
^^^^^^^^^^^^^^^^^^
MPIMod_FOUND

Imported library
^^^^^^^^^^^^^^^^
MPIMod::MPI

#]==]

set(mpi_mod_find_flags)
if(MPIMod_FIND_REQUIRED)
    list(APPEND mpi_mod_find_flags "REQUIRED")
endif()
if(MPIMod_FIND_QUIETLY)
    list(APPEND mpi_mod_find_flags "QUIET")
endif()

find_package(MPI ${mpi_mod_find_flags})

set(MPIMod_FOUND ${MPI_CXX_FOUND} CACHE BOOL "Whether or not MPI Module is found?")

if(MPIMod_FOUND AND NOT TARGET MPIMod::MPI)
    message(STATUS "Imported library MPIMod::MPI")
    message("   Include dirs: ${MPI_CXX_INCLUDE_PATH}")
    message("   Libraries: ${MPI_CXX_LIBRARIES}")
    message("   Compile flags: ${MPI_CXX_COMPILE_FLAGS}")

    if(TARGET MPI::MPI_CXX)
        add_library(MPIMod::MPI INTERFACE IMPORTED)
        target_link_libraries(MPIMod::MPI
           INTERFACE MPI::MPI_CXX
        )
    else()
        add_library(MPIMod::MPI INTERFACE IMPORTED)
        target_include_directories(MPIMod::MPI
            INTERFACE "${MPI_CXX_INCLUDE_PATH}"
        )
        target_link_libraries(MPIMod::MPI
            INTERFACE "${MPI_CXX_LIBRARIES}"
        )
    endif()

elseif(DEFINED CUSTOMIZED_MPI_DIR AND NOT TARGET MPIMod::MPI)
    message(WARNING "Unfound MPI_CXX, Please configure correct CUSTOMIZED_MPI_DIR")

    find_path(MPIMod_INCLUDE_PATH NAMES "include/mpi.h"
        PATHS ${CUSTOMIZED_MPI_DIR}
        PATH_SUFFIXES
        intel64
    )

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(MPIMod
        DEFAULT_MSG
        MPIMod_INCLUDE_PATH
    )

    add_library(MPIMod::MPI INTERFACE IMPORTED)
    target_include_directories(MPIMod::MPI
        INTERFACE "${MPIMod_INCLUDE_PATH}"
    )
    target_link_libraries(MPIMod::MPI
        INTERFACE
        libmpicxx.so libmpi.so librt.so libpthread.so libdl.so
    )
else()
    message(FATAL_ERROR "Please configure CUSTOMIZED_MPI_DIR=your_customized_path_of_mpi")
endif()
