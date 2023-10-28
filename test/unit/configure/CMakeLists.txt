################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME solver_config SOURCES test_solver_config.cpp)

################################################################################
# Tests

add_subdirectory(process)
