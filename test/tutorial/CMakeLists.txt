################################################################################
# Tests

create_standard_test(NAME README_example SOURCES test_README_example.cpp)

create_standard_test(NAME solve_results SOURCES test_solver_results.cpp)
create_standard_test(NAME multiple_grid_cells SOURCES test_multiple_grid_cells.cpp)

# These two tutorials share a mechanism whose Troe k0_A_ is a third-order rate constant converted
# from cm6 molecule-2 s-1 to MICM's mole-based units: 1.2e4 * (1e-6 * AVOGADRO)^2 = 4.35e39, which
# overflows float (FLT_MAX = 3.40e38). The parameter becomes +inf, the falloff kernel returns NaN,
# and every solve reports NaNDetected. The magnitude is inherent to the unit conversion for a
# third-order reaction, so the mechanism cannot be expressed in single precision without rescaling
# the tutorial -- and rescaling would invalidate the concentration table published alongside it in
# docs/source/user_guide/. Restrict them to double-precision builds instead.
if(MICM_USE_DOUBLE)
  create_standard_test(NAME rate_constants_no_user_defined_example_by_hand SOURCES test_rate_constants_no_user_defined_by_hand.cpp)
  create_standard_test(NAME rate_constants_user_defined_example_by_hand SOURCES test_rate_constants_user_defined_by_hand.cpp)
endif()

create_standard_test(NAME solver_configuration SOURCES test_solver_configuration.cpp)
create_standard_test(NAME vectorized_matrix_solver SOURCES test_vectorized_matrix_solver.cpp)

if(MICM_ENABLE_OPENMP)
  create_standard_test(NAME openmp_tutorial SOURCES test_openmp.cpp)
endif()

################################################################################
# Copy test data

add_custom_target(copy_tutorial_configs ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/configs ${CMAKE_BINARY_DIR}/configs)
