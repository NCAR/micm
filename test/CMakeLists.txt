# disable clang tidy for test files
set(CMAKE_CXX_CLANG_TIDY "")

add_subdirectory(unit)

if(ENABLE_REGRESSION_TESTS)
  add_subdirectory(regression)
endif()

add_subdirectory(integration)
