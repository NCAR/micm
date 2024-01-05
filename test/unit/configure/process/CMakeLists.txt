################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME arrhenius_config SOURCES test_arrhenius_config.cpp)
create_standard_test(NAME branched_config SOURCES test_branched_config.cpp)
create_standard_test(NAME emission_config SOURCES test_emission_config.cpp)
create_standard_test(NAME first_order_loss_config SOURCES test_first_order_loss_config.cpp)
create_standard_test(NAME photolysis_config SOURCES test_photolysis_config.cpp)
create_standard_test(NAME surface_config SOURCES test_surface_config.cpp)
create_standard_test(NAME ternary_chemical_activation_config SOURCES test_ternary_chemical_activation_config.cpp)
create_standard_test(NAME troe_config SOURCES test_troe_config.cpp)
create_standard_test(NAME tunneling_config SOURCES test_tunneling_config.cpp)
create_standard_test(NAME user_defined_config SOURCES test_user_defined_config.cpp)
