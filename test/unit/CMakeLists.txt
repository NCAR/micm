if(NOT ${MICM_GPU_TYPE} STREQUAL "None")
  add_subdirectory(cuda)
endif()
if(MICM_ENABLE_LLVM)
  add_subdirectory(jit)
endif()
add_subdirectory(process)
add_subdirectory(solver)
add_subdirectory(system)
add_subdirectory(util)
add_subdirectory(version)
if(MICM_ENABLE_OPENMP)
  add_subdirectory(openmp)
endif()