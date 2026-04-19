################################################################################
# Tests

create_standard_test(NAME chapman_integration SOURCES test_chapman_integration.cpp)
create_standard_test(NAME analytical_rosenbrock_integration SOURCES test_analytical_rosenbrock.cpp)
create_standard_test(NAME analytical_backward_euler SOURCES test_analytical_backward_euler.cpp)
create_standard_test(NAME terminator SOURCES test_terminator.cpp)
create_standard_test(NAME integrated_reaction_rates SOURCES test_integrated_reaction_rates.cpp)
create_standard_test(NAME aerosol_model_backward_euler SOURCES test_aerosol_model_backward_euler.cpp)
create_standard_test(NAME aerosol_model_rosenbrock SOURCES test_aerosol_model_rosenbrock.cpp)
create_standard_test(NAME equilibrium_constraint SOURCES test_equilibrium_constraint.cpp)
create_standard_test(NAME linear_constraint SOURCES test_linear_constraint.cpp)
create_standard_test(NAME external_model_constraints SOURCES test_external_model_constraints.cpp)
create_standard_test(NAME dae_constraint_overshoot SOURCES test_dae_constraint_overshoot.cpp)

if(NOT ${MICM_GPU_TYPE} STREQUAL "None")
  add_subdirectory(cuda)
endif()
