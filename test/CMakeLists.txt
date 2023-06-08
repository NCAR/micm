# disable clang tidy for test files
set(CMAKE_CXX_CLANG_TIDY "")

if(ENABLE_LOCAL_TESTS)
  add_subdirectory(unit)

  if(ENABLE_REGESSION_TESTS)
    add_subdirectory(regression)
  endif()

  add_subdirectory(integration)
endif()

if(ENABLE_HPC_TESTS)
  add_subdirectory(hpc_tests)
endif()
