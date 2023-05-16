################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME arrhenius_rate_constant SOURCES test_arrhenius_rate_constant.cpp)
create_standard_test(NAME photolysis_rate_constant SOURCES test_photolysis_rate_constant.cpp)
create_standard_test(NAME troe_rate_constant SOURCES test_troe_rate_constant.cpp)
