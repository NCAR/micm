################################################################################
# Copy test data

add_custom_target(copy_example_configs ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_SOURCE_DIR}/examples ${CMAKE_BINARY_DIR}/examples)

################################################################################
# Build example tests

create_standard_test(NAME Chapman SOURCES test_Chapman.cpp)
create_standard_test(NAME TS1 SOURCES test_TS1.cpp)
