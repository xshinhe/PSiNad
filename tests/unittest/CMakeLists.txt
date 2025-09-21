include(CTest)
include(Catch)

# Find Catch2 library (prefer system-installed version)
find_package(Catch2 3 QUIET)

# Add include directories (for bundled Catch2)
include_directories(
    ${CMAKE_SOURCE_DIR}/thirdpart/catch2/src
)

# Get all test source files
file(GLOB test_list *.cpp)

# Set C++ standard (keep consistent with main project)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## Build each test target
foreach(test_file ${test_list})
    # Extract filename without path and extension
    get_filename_component(file_name ${test_file} NAME_WE)
    set(test_name ${file_name})

    # Create test executable
    add_executable(${test_name} ${test_file})

    # Link necessary libraries - Include Catch2WithMain for auto-generated main()
    target_link_libraries(${test_name} 
        PRIVATE 
        Catch2::Catch2 
        Catch2::Catch2WithMain  # Provides required main() function
        ${SHARED_TARGET}
    )

    # Register test with CTest
    catch_discover_tests(${test_name})

    # Additional linking settings for macOS arm64 architecture
    if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
        target_link_options(${test_name} PRIVATE "-arch" "arm64")
    endif()
endforeach()
