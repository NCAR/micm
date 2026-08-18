################################################################################
# Tests

create_standard_test(NAME kokkos_constraint_set SOURCES test_kokkos_constraint_set.cpp LIBRARIES musica::micm_kokkos)

add_subdirectory(type)