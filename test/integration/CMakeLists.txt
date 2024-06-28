################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME chapman_integration SOURCES chapman.cpp)
create_standard_test(NAME analytical_rosenbrock_integration SOURCES analytical_rosenbrock.cpp)
create_standard_test(NAME terminator_integration SOURCES terminator.cpp)
create_standard_test(NAME integrated_reaction_rates SOURCES integrated_reaction_rates.cpp)

if(MICM_ENABLE_LLVM)
  add_subdirectory(jit)
endif()

if(MICM_ENABLE_CUDA)
  add_subdirectory(cuda)
endif()
