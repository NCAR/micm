################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

add_library(micm_cuda_test_utils)
target_link_libraries(micm_cuda_test_utils PUBLIC cudart micm)
target_sources(micm_cuda_test_utils PRIVATE cuda_matrix_utils.cu)
set_target_properties(micm_cuda_test_utils PROPERTIES LINKER_LANGUAGE CXX)

# add a library to use customized GTest main function for CUDA tests
add_library(cuda_gtest_main STATIC cuda_gtest_main.cpp)
target_link_libraries(cuda_gtest_main PUBLIC gtest micm)

create_standard_test(NAME cuda_dense_matrix SOURCES test_cuda_dense_matrix.cpp LIBRARIES musica::micm_cuda micm_cuda_test_utils IS_CUDA_TEST)
create_standard_test(NAME cuda_sparse_matrix SOURCES test_cuda_sparse_matrix.cpp LIBRARIES musica::micm_cuda micm_cuda_test_utils IS_CUDA_TEST)