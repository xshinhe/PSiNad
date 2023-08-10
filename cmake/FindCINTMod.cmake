#[==[

FindCINTMod
-----------

Imported library
^^^^^^^^^^^^^^^^
CINTMod::CINT

#]==]


find_library(LIBCINT libcint.so ${LD_LIBRARY_PATH})


if(LIBCINT)
    option(USE_LIBCINT "use libcint library" ON)

    add_library(CINTMod::CINT INTERFACE IMPORTED)

    target_link_libraries(CINTMod::CINT
        INTERFACE "${LIBCINT}"
    )
else()
    option(USE_LIBCINT "use libcint library" OFF)
    add_library(CINTMod::CINT INTERFACE IMPORTED)

    target_link_libraries(CINTMod::CINT
        INTERFACE ""
    )
endif(LIBCINT)