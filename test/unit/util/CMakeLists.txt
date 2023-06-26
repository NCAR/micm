################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME matrix SOURCES test_matrix.cpp)
create_standard_test(NAME sparse_matrix SOURCES test_sparse_matrix.cpp)
create_standard_test(NAME vector_matrix SOURCES test_vector_matrix.cpp)
