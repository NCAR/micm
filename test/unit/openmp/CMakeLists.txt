################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME openmp SOURCES test_openmp.cpp)