################################################################################
# Tests

create_standard_test(NAME kokkos_equilibrium SOURCES test_kokkos_equilibrium.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_linear SOURCES test_kokkos_linear.cpp LIBRARIES musica::micm_kokkos)