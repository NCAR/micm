################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME chapman_integration SOURCES chapman.cpp)
create_standard_test(NAME analytical_rosenbrock_integration SOURCES analytical_rosenbrock.cpp)
create_standard_test(NAME terminator_integration SOURCES terminator.cpp)

if(MICM_ENABLE_LLVM)
  create_standard_test(NAME analytical_jit_rosenbrock_integration SOURCES analytical_jit_rosenbrock.cpp)
  create_standard_test(NAME jit_terminator_integration SOURCES jit_terminator.cpp)
endif()
