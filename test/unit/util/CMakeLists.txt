################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME matrix SOURCES test_matrix.cpp)
create_standard_test(NAME sparse_matrix_standard_ordering SOURCES test_sparse_matrix_standard_ordering.cpp)
create_standard_test(NAME sparse_matrix_vector_ordering SOURCES test_sparse_matrix_vector_ordering.cpp)
create_standard_test(NAME vector_matrix SOURCES test_vector_matrix.cpp)

add_library(micm_cuda_test_utils)
target_link_libraries(micm_cuda_test_utils PUBLIC cudart micm)
target_sources(micm_cuda_test_utils PRIVATE cuda_matrix_utils.cu)

if(ENABLE_CUDA)
  create_standard_test(NAME cuda_vector_matrix SOURCES test_cuda_vector_matrix.cpp LIBRARIES musica::micm_cuda micm_cuda_test_utils)
endif()
