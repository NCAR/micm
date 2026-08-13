################################################################################
# Tests

add_executable(test_kokkos_analytical_backward_euler test_kokkos_analytical_backward_euler.cpp)
target_link_libraries(test_kokkos_analytical_backward_euler PUBLIC musica::micm musica::micm_kokkos GTest::gtest)
add_test(NAME kokkos_analytical_backward_euler COMMAND test_kokkos_analytical_backward_euler WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_kokkos_analytical_rosenbrock test_kokkos_analytical_rosenbrock.cpp)
target_link_libraries(test_kokkos_analytical_rosenbrock PUBLIC musica::micm musica::micm_kokkos GTest::gtest)
add_test(NAME kokkos_analytical_rosenbrock COMMAND test_kokkos_analytical_rosenbrock WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_kokkos_chapman_integration test_kokkos_chapman_integration.cpp)
target_link_libraries(test_kokkos_chapman_integration PUBLIC musica::micm musica::micm_kokkos GTest::gtest)
add_test(NAME kokkos_chapman_integration COMMAND test_kokkos_chapman_integration WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_kokkos_dae_algebraic_error_insensitivity test_kokkos_dae_algebraic_error_insensitivity.cpp)
target_link_libraries(test_kokkos_dae_algebraic_error_insensitivity PUBLIC musica::micm musica::micm_kokkos GTest::gtest)
add_test(NAME kokkos_dae_algebraic_error_insensitivity COMMAND test_kokkos_dae_algebraic_error_insensitivity WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_kokkos_dae_constraint_overshoot test_kokkos_dae_constraint_overshoot.cpp)
target_link_libraries(test_kokkos_dae_constraint_overshoot PUBLIC musica::micm musica::micm_kokkos GTest::gtest)
add_test(NAME kokkos_dae_constraint_overshoot COMMAND test_kokkos_dae_constraint_overshoot WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_kokkos_equilibrium_constraint test_kokkos_equilibrium_constraint.cpp)
target_link_libraries(test_kokkos_equilibrium_constraint PUBLIC musica::micm musica::micm_kokkos GTest::gtest)
add_test(NAME kokkos_equilibrium_constraint COMMAND test_kokkos_equilibrium_constraint WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_kokkos_linear_constraint test_kokkos_linear_constraint.cpp)
target_link_libraries(test_kokkos_linear_constraint PUBLIC musica::micm musica::micm_kokkos GTest::gtest)
add_test(NAME kokkos_linear_constraint COMMAND test_kokkos_linear_constraint WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_kokkos_terminator test_kokkos_terminator.cpp)
target_link_libraries(test_kokkos_terminator PUBLIC musica::micm musica::micm_kokkos GTest::gtest)
add_test(NAME kokkos_terminator COMMAND test_kokkos_terminator WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
