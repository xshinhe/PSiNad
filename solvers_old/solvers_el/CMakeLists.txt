file(GLOB_RECURSE el_srcs *.cpp *.c)

if(DEFINED USE_LIBCINT)
    if(USE_LIBCINT)
        add_definitions(-DUSE_LIBCINT)
    endif()
endif()

add_library(el ${el_srcs})

target_link_libraries(
    el
    utils
    base
    models
    glog::glog
    CINTMod::CINT
)

target_precompile_headers(el REUSE_FROM utils)

