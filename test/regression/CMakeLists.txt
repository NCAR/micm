################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME regression SOURCES test_regression.cpp)

add_subdirectory(RosenbrockChapman)

################################################################################
# Copy expected results

add_custom_target(copy_expected_regression_results ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/expected_results ${CMAKE_BINARY_DIR}/expected_results)