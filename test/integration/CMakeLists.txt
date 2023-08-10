################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME chapman_integration SOURCES chapman.cpp)
create_standard_test(NAME analytical_integration SOURCES analytical.cpp)