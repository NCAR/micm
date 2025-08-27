################################################################################
# Tests

create_standard_test(NAME chapman_integration SOURCES test_chapman_integration.cpp)
create_standard_test(NAME analytical_rosenbrock_integration SOURCES test_analytical_rosenbrock.cpp)
create_standard_test(NAME analytical_backward_euler SOURCES test_analytical_backward_euler.cpp)
create_standard_test(NAME terminator SOURCES test_terminator.cpp)
create_standard_test(NAME integrated_reaction_rates SOURCES test_integrated_reaction_rates.cpp)

if(MICM_ENABLE_LLVM)
  add_subdirectory(jit)
endif()

if(NOT ${MICM_GPU_TYPE} STREQUAL "None")
  add_subdirectory(cuda)
endif()
