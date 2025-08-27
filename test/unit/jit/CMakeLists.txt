################################################################################
# Tests

create_standard_test(NAME jit_function SOURCES test_jit_function.cpp)
add_subdirectory(process)
add_subdirectory(solver)