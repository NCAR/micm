################################################################################
# Tests

create_standard_test(NAME README_example SOURCES test_README_example.cpp)

create_standard_test(NAME solve_results SOURCES test_solver_results.cpp)
create_standard_test(NAME multiple_grid_cells SOURCES test_multiple_grid_cells.cpp)
create_standard_test(NAME rate_constants_no_user_defined_example_by_hand SOURCES test_rate_constants_no_user_defined_by_hand.cpp)
create_standard_test(NAME rate_constants_user_defined_example_by_hand SOURCES test_rate_constants_user_defined_by_hand.cpp)
create_standard_test(NAME solver_configuration SOURCES test_solver_configuration.cpp)
create_standard_test(NAME vectorized_matrix_solver SOURCES test_vectorized_matrix_solver.cpp)

if(MICM_ENABLE_OPENMP)
  create_standard_test(NAME openmp_tutorial SOURCES test_openmp.cpp)
endif()

################################################################################
# Copy test data

add_custom_target(copy_tutorial_configs ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/configs ${CMAKE_BINARY_DIR}/configs)
