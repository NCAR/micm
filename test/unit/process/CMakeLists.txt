################################################################################
# Tests
create_standard_test(NAME process SOURCES test_process.cpp)
create_standard_test(NAME process_set SOURCES test_process_set.cpp)
create_standard_test(NAME reaction_rate_store SOURCES test_reaction_rate_store.cpp)

add_subdirectory(rate_constant)
add_subdirectory(transfer_coefficient)