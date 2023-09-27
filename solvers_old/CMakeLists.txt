add_library(base solver.cpp)
target_link_libraries(
    base
    utils
    models
    glog::glog
    MPIMod::MPI
    HighFive
)

# configure solvers needed for compilation
if(EXISTS ${CMAKE_BINARY_DIR}/config.json)
    set(config_file ${CMAKE_BINARY_DIR}/config.json) # customized config.json!
else()
    set(config_file ${CMAKE_SOURCE_DIR}/config.json) # default with all solvers
endif()

execute_process(
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/generate_factory.py
            ${CMAKE_SOURCE_DIR}
            ${config_file}
            solvers
            ${CMAKE_SOURCE_DIR}/solvers/solverfactory.cpp
    TIMEOUT 5
    RESULT_VARIABLE _status
    OUTPUT_VARIABLE specify_solvers_srcs
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(test_list ${specify_solvers_srcs})
list (LENGTH test_list leng1)
foreach(sflag md PI nad rd el dev)

    list (FILTER test_list EXCLUDE REGEX "solvers_${sflag}")
    list (LENGTH test_list leng2)

    if(${leng2} LESS ${leng1})
        message(STATUS "adding solvers_${sflag}")
        add_subdirectory(solvers_${sflag})
        set(solvers_added ${solvers_added} ${sflag})
    endif()
    set(leng1 ${leng2})
endforeach()

add_library(solvers solverfactory.cpp)
target_link_libraries(
    solvers
    utils
    models
    base
    ${solvers_added}
    glog::glog
    MPIMod::MPI
    HighFive
)

target_precompile_headers(solvers REUSE_FROM utils)
