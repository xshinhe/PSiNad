# python
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -m pybind11 --includes
    TIMEOUT 5
    RESULT_VARIABLE _status
    OUTPUT_VARIABLE pybind11_include_flag
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(INFO "  using pybind11 flags: ${pybind11_include_flag}")
add_definitions("${pybind11_include_flag}")

# configure source file
# python
if(EXISTS ${CMAKE_BINARY_DIR}/config.json)
    set(config_file ${CMAKE_BINARY_DIR}/config.json)
else()
    set(config_file ${CMAKE_SOURCE_DIR}/config.json)
endif()

execute_process(
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/generate_factory.py
            ${CMAKE_SOURCE_DIR}
            ${config_file}
            models
            ${CMAKE_SOURCE_DIR}/models/modelfactory.cpp
    TIMEOUT 5
    RESULT_VARIABLE _status
    OUTPUT_VARIABLE specify_models_srcs_raw
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

file(GLOB_RECURSE models_srcs
    # intrinsic
    bath/*.cpp
    thermostat/*.cpp
    model.cpp
    forcefieldbase.cpp
    modelfactory.cpp
)
string(REPLACE "\n" ";"  specify_models_srcs ${specify_models_srcs_raw})
list(APPEND models_srcs ${specify_models_srcs})
add_library(models ${models_srcs})


set(test_copy ${specify_models_srcs})
list (FILTER test_copy INCLUDE REGEX py_interface)
list (LENGTH test_copy leng)
if(${leng} EQUAL 0)
    set(whetherpybind "")
else()
    set(whetherpybind pybind11::embed)
endif()

set(test_copy ${specify_models_srcs})
list (FILTER test_copy INCLUDE REGEX bo_forcefield)
list (LENGTH test_copy leng)
if(${leng} EQUAL 0)
    set(whetherflib ${whetherflib})
else()
    # fortran code library
    add_subdirectory(flib)
    set(whetherflib h2o h2o2 nh3 parah2 qspcfw liquidne)
endif()

set(test_copy ${specify_models_srcs})
list (FILTER test_copy INCLUDE REGEX sanderapi)
list (LENGTH test_copy leng)
if(${leng} EQUAL 0)
else()
    ### blame ###
    # DON'T WRITE YOUR LOCAL PATH INTO NON-ROOT-DIR CMAKELISTS.TXT; pass it by cmake -DAMBERHOME=xxxxx
    # set(AMBERHOME /home/zhf/work/amberfold/amber20/)
    set(whetherflib ${whetherflib} sander)
    include_directories(${AMBERHOME}/include/)
    link_directories(${AMBERHOME}/lib/)
endif()



target_link_libraries(
    models
    utils
    glog::glog
    # other libraries
    ${whetherflib}
    ${whetherpybind}
)

target_precompile_headers(models REUSE_FROM utils)
