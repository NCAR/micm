if(MICM_ENABLE_JSON)
  add_subdirectory(configure)
endif()
if(MICM_ENABLE_LLVM)
  add_subdirectory(jit)
endif()
add_subdirectory(process)
add_subdirectory(solver)
add_subdirectory(system)
add_subdirectory(util)
add_subdirectory(version)
if(MICM_ENABLE_OPENMP AND MICM_ENABLE_JSON)
  add_subdirectory(openmp)
endif()

################################################################################
# Copy test data

add_custom_target(copy_unit_test_configs ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/unit_configs ${CMAKE_BINARY_DIR}/unit_configs)