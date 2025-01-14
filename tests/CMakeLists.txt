include(CTest)
include(Catch)

file(GLOB test_list unittest/*.cpp)

## Base tests
foreach(test_file ${test_list})
  string(REGEX MATCH "[A-Za-z0-9_]*.cpp" file_name ${test_file})
  string(REPLACE ".cpp" "" test_name ${file_name})
  add_executable(${test_name} ${file_name})
  target_link_libraries(${test_name} PRIVATE Catch2::Catch2)
  catch_discover_tests(${test_name})
endforeach()
