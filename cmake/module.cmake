macro(module_config)
	find_package(PythonInterp REQUIRED) # help for configuration

	# configure solvers needed for compilation
	if(EXISTS ${CMAKE_BINARY_DIR}/config.json)
	    set(config_file ${CMAKE_BINARY_DIR}/config.json) # customized config.json!
	else()
	    set(config_file ${CMAKE_SOURCE_DIR}/config.json) # default with all solvers
	endif()

	execute_process(
	    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/generate_helper/generate_factory.py
	            ${CMAKE_SOURCE_DIR}
	            ${config_file}
	            solvers
	            ${CMAKE_SOURCE_DIR}/generate/solverfactory.cpp
	    TIMEOUT 5
	    RESULT_VARIABLE _status
	    OUTPUT_VARIABLE specify_solvers_srcs
	    OUTPUT_STRIP_TRAILING_WHITESPACE
	)

	execute_process(
	    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/generate_helper/generate_factory.py
	            ${CMAKE_SOURCE_DIR}
	            ${config_file}
	            models
	            ${CMAKE_SOURCE_DIR}/generate/modelfactory.cpp
	    TIMEOUT 5
	    RESULT_VARIABLE _status
	    OUTPUT_VARIABLE specify_models_srcs_raw
	    OUTPUT_STRIP_TRAILING_WHITESPACE
	)

endmacro()