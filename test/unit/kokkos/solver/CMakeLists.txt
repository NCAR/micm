################################################################################
# Tests

create_standard_test(NAME kokkos_lu_decomposition_mozart_in_place SOURCES test_kokkos_lu_decomposition_mozart_in_place.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_linaer_solver_in_place SOURCES test_kokkos_linear_solver_in_place.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_backward_euler SOURCES test_kokkos_backward_euler.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_rosenbrock SOURCES test_kokkos_rosenbrock.cpp LIBRARIES musica::micm_kokkos)
