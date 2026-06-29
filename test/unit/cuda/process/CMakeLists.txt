################################################################################
# Tests

create_standard_test(NAME cuda_process_set SOURCES test_cuda_process_set.cpp LIBRARIES musica::micm_cuda IS_CUDA_TEST)
create_standard_test(NAME cuda_rate_constants SOURCES test_cuda_rate_constants.cpp LIBRARIES musica::micm_cuda IS_CUDA_TEST)