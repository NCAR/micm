################################################################################
# Tests

create_standard_test(NAME cuda_process_set SOURCES test_cuda_process_set.cpp LIBRARIES musica::micm_cuda IS_CUDA_TEST)
create_standard_test(NAME cuda_SubtractJacobianTerms SOURCES test_cuda_SubtractJacobianTerms.cpp LIBRARIES musica::micm_cuda IS_CUDA_TEST)