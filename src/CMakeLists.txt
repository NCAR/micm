################################################################################
# micm version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/micm/version.hpp @ONLY)

################################################################################
# micm solver
add_library(micm INTERFACE)
add_library(musica::micm ALIAS micm)

# target_compile_features(micm INTERFACE cxx_std_20)

target_include_directories(micm
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(ENABLE_JSON)
  target_link_libraries(micm INTERFACE nlohmann_json::nlohmann_json)
  target_compile_definitions(micm INTERFACE USE_JSON)
endif()

if(ENABLE_OPENMP)
  target_link_libraries(micm INTERFACE OpenMP::OpenMP_CXX)
endif()

if(ENABLE_MPI)
  target_link_libraries(micm INTERFACE MPI::MPI_CXX)
endif()

if(ENABLE_CUDA OR ENABLE_OPENACC)
  add_library(micm_gpu)
  add_library(musica::micm_gpu ALIAS micm_gpu)

  if(NOT GPU_TYPE)
    message(FATAL_ERROR "GPU_TYPE is not set or is empty. Please provide a GPU type.")
  endif()

  set(GPU_ARCH "")
  if (GPU_TYPE STREQUAL a100)
    set(GPU_ARCH "cc80")
  endif()
  if (GPU_TYPE STREQUAL v100)
    set(GPU_ARCH "cc70")
  endif()

  if(NOT GPU_ARCH)
    message(FATAL_ERROR "GPU_TYPE (${GPU_TYPE}) is not recognized. Available options are [a100, v100].")
  endif()
  
  set(GPU_FLAGS -gpu=${GPU_ARCH},lineinfo,nofma -Minfo=accel)
  message(STATUS "GPU Flags: ${OpenACC_CXX_FLAGS}  ${GPU_FLAGS}")

  target_compile_options(micm_gpu PRIVATE ${GPU_FLAGS})

  target_link_libraries(micm_gpu 
    PRIVATE micm
  )

  if(ENABLE_CUDA)
    target_link_libraries(micm_gpu 
      PUBLIC cudart
    )
  endif()

  if(ENABLE_OPENACC)
    target_link_libraries(micm_gpu 
      PUBLIC OpenACC::OpenACC_CXX
    )
  endif()

  add_subdirectory(process)
endif()
