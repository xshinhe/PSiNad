#[==[

FindMKLMod
-----------

Resulted variables
^^^^^^^^^^^^^^^^^^
MKLMod_FOUND

Imported library
^^^^^^^^^^^^^^^^
MKLMod::MKL

#]==]

set(mkl_mod_find_flags)
if(MKLMod_FIND_REQUIRED)
    list(APPEND mkl_mod_find_flags "REQUIRED")
else()
    list(APPEND mkl_mod_find_flags "QUIET")
endif()

set(MKL_INTERFACE "lp64")
set(MKL_LINK "dynamic")
find_package(MKL ${mkl_mod_find_flags})

set(MKLMod_FOUND ${MKL_FOUND} CACHE BOOL "Whether or not MKL Module is found?")

if(MKLMod_FOUND AND NOT TARGET MKLMod::MKL)
    message(STATUS "Imported library MKLMod::MKL")
    message("   Include dirs: ${MKL_INCLUDE_PATH}")
    message("   Libraries: ${MKL_LIBRARIES}")
    message("   Comklle flags: ${MKL_COMPILE_FLAGS}")

    if(TARGET MKL::MKL)
        add_library(MKLMod::MKL INTERFACE IMPORTED)
        target_link_libraries(MKLMod::MKL
            INTERFACE MKL::MKL
        )
    else()
        add_library(MKLMod::MKL INTERFACE IMPORTED)
        target_include_directories(MKLMod::MKL
            INTERFACE "${MKL_INCLUDE_PATH}"
        )
        target_link_libraries(MKLMod::MKL
            INTERFACE "${MKL_LIBRARIES}"
        )
    endif()

elseif(CUSTOMIZED_MKL_DIR AND NOT TARGET MKLMod::MKL)
    find_path(MKLMod_INCLUDE_PATH NAMES "mkl.h"
        PATHS ${CUSTOMIZED_MKL_DIR}
        PATH_SUFFIXES
        include
    )

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(MKLMod
        DEFAULT_MSG
        MKLMod_INCLUDE_PATH
    )

    set(MKLMod_LIBRARIES
        libiomp5.so
        libmkl_core.so
        libmkl_rt.so
        libmkl_intel_thread.so
        libmkl_intel_lp64.so
    )

    add_library(MKLMod::MKL INTERFACE IMPORTED)

    target_include_directories(MKLMod::MKL
        INTERFACE "${MKLMod_INCLUDE_PATH}"
    )
    target_link_directories(MKLMod::MKL
        INTERFACE "${CUSTOMIZED_MKL_DIR}"
        INTERFACE "${CUSTOMIZED_MKL_DIR}/lib"
        INTERFACE "${CUSTOMIZED_MKL_DIR}/lib/intel64"
    )
    target_link_libraries(MKLMod::MKL
        INTERFACE "${MKLMod_LIBRARIES}"
    )


else()
    message(FATAL_ERROR "Please configure CUSTOMIZED_MKL_DIR=your_customized_path_of_mkl")
endif()