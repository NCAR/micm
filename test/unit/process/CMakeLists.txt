################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME process SOURCES test_process.cpp)
create_standard_test(NAME process_set SOURCES test_process_set.cpp)
create_standard_test(NAME arrhenius_rate_constant SOURCES test_arrhenius_rate_constant.cpp)
create_standard_test(NAME branched_rate_constant SOURCES test_branched_rate_constant.cpp)
create_standard_test(NAME lambda_rate_constant SOURCES test_lambda_rate_constant.cpp)
create_standard_test(NAME surface_rate_constant SOURCES test_surface_rate_constant.cpp)
create_standard_test(NAME taylor_series_rate_constant SOURCES test_taylor_series_rate_constant.cpp)
create_standard_test(NAME ternary_chemical_activation_rate_constant SOURCES test_ternary_chemical_activation_rate_constant.cpp)
create_standard_test(NAME troe_rate_constant SOURCES test_troe_rate_constant.cpp)
create_standard_test(NAME tunneling_rate_constant SOURCES test_tunneling_rate_constant.cpp)
create_standard_test(NAME user_defined_rate_constant SOURCES test_user_defined_rate_constant.cpp)