# disable clang tidy for test files
set(CMAKE_CXX_CLANG_TIDY "")

################################################################################
# Test utilities

include(test_util)

################################################################################
if(MICM_ENABLE_BENCHMARK)
  add_subdirectory(benchmark)
endif()
add_subdirectory(unit)
add_subdirectory(integration)
add_subdirectory(regression)
add_subdirectory(tutorial)
