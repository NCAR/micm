# disable clang tidy for test files
set(CMAKE_CXX_CLANG_TIDY "")

if(ENABLE_LOCAL_TESTS)
  add_subdirectory(unit)

  if(ENABLE_REGESSION_TESTS)
    add_subdirectory(regression)
  endif()

  add_subdirectory(integration)
endif()
