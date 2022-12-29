
################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME arrhenius_rate_constant SOURCES test_arrhenius_rate_constant.cpp)
create_standard_test(NAME external_rate_constant SOURCES test_external_rate_constant.cpp)
