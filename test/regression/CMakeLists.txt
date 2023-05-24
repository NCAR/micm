################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

# compile the old micm library
add_subdirectory(old_micm)
add_subdirectory(HardCodedChapman)
add_subdirectory(RosenbrockChapman)

################################################################################
# Copy test data

add_custom_target(copy_regression_test_configs ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/regression_configs ${CMAKE_BINARY_DIR}/regression_configs)