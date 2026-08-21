################################################################################
# Tests

create_standard_test(NAME kokkos_analytical_backward_euler SOURCES test_kokkos_analytical_backward_euler.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_analytical_rosenbrock SOURCES test_kokkos_analytical_rosenbrock.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_chapman_integration SOURCES test_kokkos_chapman_integration.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_cpu_agreement SOURCES test_kokkos_cpu_agreement.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_dae_algebraic_error_insensitivity SOURCES test_kokkos_dae_algebraic_error_insensitivity.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_dae_constraint_overshoot SOURCES test_kokkos_dae_constraint_overshoot.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_equilibrium_constraint SOURCES test_kokkos_equilibrium_constraint.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_linear_constraint SOURCES test_kokkos_linear_constraint.cpp LIBRARIES musica::micm_kokkos)
create_standard_test(NAME kokkos_terminator SOURCES test_kokkos_terminator.cpp LIBRARIES musica::micm_kokkos)