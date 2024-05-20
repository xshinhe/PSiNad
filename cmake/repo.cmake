# get git hash
macro(get_repo_info _repo_branch _repo_date _repo_hash _repo_status)
    find_package(Git QUIET)
    if(GIT_FOUND)
      execute_process(
        COMMAND ${GIT_EXECUTABLE} status -s
        OUTPUT_VARIABLE ${_repo_status}
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        WORKING_DIRECTORY
          ${CMAKE_SOURCE_DIR}
        )

      execute_process(
        COMMAND ${GIT_EXECUTABLE} symbolic-ref --short -q HEAD
        OUTPUT_VARIABLE ${_repo_branch}
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        WORKING_DIRECTORY
          ${CMAKE_SOURCE_DIR}
        )

      execute_process(
        COMMAND ${GIT_EXECUTABLE} log -1 --format=%cd --date=short
        OUTPUT_VARIABLE ${_repo_date}
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        WORKING_DIRECTORY
          ${CMAKE_SOURCE_DIR}
        )

      execute_process(
        COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%H
        OUTPUT_VARIABLE ${_repo_hash}
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        WORKING_DIRECTORY
          ${CMAKE_SOURCE_DIR}
        )
    endif()
    if(NOT ${_repo_hash})
      set(${_repo_hash} "Unknown")
    endif()
endmacro()
