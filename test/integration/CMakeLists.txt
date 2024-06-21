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
  create_standard_test(NAME analytical_jit_rosenbrock_integration SOURCES analytical_jit_rosenbrock.cpp)
  create_standard_test(NAME jit_terminator_integration SOURCES jit_terminator.cpp)
endif()

if(MICM_ENABLE_CUDA)
  # the compute-sanitizer hangs for any integration test somehow; turns it off for now
  create_standard_test(NAME analytical_cuda_rosenbrock_integration SOURCES analytical_cuda_rosenbrock.cpp LIBRARIES musica::micm_cuda)
endif()
