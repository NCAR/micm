################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME regression_test_rosenbrock_p_force SOURCES regression_test_p_force.cpp)
create_standard_test(NAME regression_test_rosenbrock_dforce_dy SOURCES regression_test_dforce_dy.cpp)