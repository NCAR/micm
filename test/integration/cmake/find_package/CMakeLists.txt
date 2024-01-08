cmake_minimum_required(VERSION 3.11)

project(
  test_micm_fetch_content
  VERSION 0.0.0
  LANGUAGES CXX
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

find_package(micm 3.3.1 REQUIRED)

################################################################################
# Tests
add_executable(test_micm ../test_micm.cpp)

target_link_libraries(test_micm 
  PUBLIC 
    musica::micm
)

enable_testing()

add_test(
  NAME test_micm 
  COMMAND test_micm
)

add_custom_target(copy_configs ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/../configs ${CMAKE_BINARY_DIR}/configs)