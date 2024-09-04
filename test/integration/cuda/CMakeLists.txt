################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME cuda_analytical_rosenbrock_integration SOURCES test_cuda_analytical_rosenbrock.cpp LIBRARIES musica::micm_cuda IS_CUDA_TEST)