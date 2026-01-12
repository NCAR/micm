################################################################################
# Tests
create_standard_test(NAME aerosol_scope SOURCES test_process_configuration.cpp)
create_standard_test(NAME process SOURCES test_process.cpp)
create_standard_test(NAME process_set SOURCES test_process_set.cpp)

add_subdirectory(rate_constant)
add_subdirectory(transfer_coefficient)