# disable clang tidy for test files
set(CMAKE_CXX_CLANG_TIDY "")

add_executable(micmDriver example.cpp)
target_link_libraries(micmDriver PUBLIC musica::micm )
################################################################################
# Copy test data

add_custom_target(copy_example_test_configs ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/configs ${CMAKE_BINARY_DIR}/examples/configs)

