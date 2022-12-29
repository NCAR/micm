
################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME arrhenius_rate_constant SOURCES test_arrhenius_rate_constant.cpp)
create_standard_test(NAME external_rate_constant SOURCES test_external_rate_constant.cpp)
create_standard_test(NAME intraphase_process SOURCES test_intraphase_process.cpp)
create_standard_test(NAME intraphase_process_builder SOURCES test_intraphase_process_builder.cpp)
