################################################################################
# Tests

create_standard_test(NAME kokkos_dense_matrix SOURCES test_kokkos_dense_matrix.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_padded_vector SOURCES test_kokkos_padded_vector.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_scalar_view SOURCES test_kokkos_scalar_view.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_sparse_matrix SOURCES test_kokkos_sparse_matrix.cpp LIBRARIES musica::micm_kokkos)
