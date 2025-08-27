################################################################################
# Tests

create_standard_test(NAME matrix SOURCES test_matrix.cpp)
create_standard_test(NAME sparse_matrix_standard_compressed_column_ordering SOURCES test_sparse_matrix_standard_ordering_column.cpp)
create_standard_test(NAME sparse_matrix_standard_compressed_row_ordering SOURCES test_sparse_matrix_standard_ordering_row.cpp)
create_standard_test(NAME sparse_matrix_vector_compressed_column_ordering SOURCES test_sparse_matrix_vector_ordering_column.cpp)
create_standard_test(NAME sparse_matrix_vector_compressed_row_ordering SOURCES test_sparse_matrix_vector_ordering_row.cpp)
create_standard_test(NAME vector_matrix SOURCES test_vector_matrix.cpp)