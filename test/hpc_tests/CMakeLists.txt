################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_hpc_test(NAME chapman_integration SOURCES chapman_integration.cpp)