################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME branched_config SOURCES test_branched_config.cpp)
create_standard_test(NAME ternary_chemical_activation_config SOURCES test_ternary_chemical_activation_config.cpp)
create_standard_test(NAME tunneling_config SOURCES test_tunneling_config.cpp)