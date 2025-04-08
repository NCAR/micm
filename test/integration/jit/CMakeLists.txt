################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME jit_analytical_rosenbrock_integration SOURCES test_jit_analytical_rosenbrock.cpp)
create_standard_test(NAME jit_terminator_integration SOURCES test_jit_terminator.cpp)