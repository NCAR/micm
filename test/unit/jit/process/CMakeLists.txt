################################################################################
# Tests

create_standard_test(NAME jit_process_set SOURCES test_jit_process_set.cpp)
create_standard_test(NAME jit_forcing_calculation SOURCES test_jit_forcing_calculation.cpp forcing_calculation.cpp)