################################################################################
# micm version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/micm/version.hpp @ONLY)

################################################################################
# micm solver
add_library(micm INTERFACE)
add_library(musica::micm ALIAS micm)

target_compile_features(micm INTERFACE cxx_std_20)

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

if(ENABLE_OPENACC OR ENABLE_CUDA)
  if(NOT GPU_TYPE)
    message(FATAL_ERROR "GPU_TYPE is not set or is empty. Please provide a GPU type.")
  endif()

  set(GPU_ARCH "")
  if (GPU_TYPE STREQUAL a100)
    set(GPU_ARCH "80")
  endif()
  if (GPU_TYPE STREQUAL v100)
    set(GPU_ARCH "70")
  endif()

  if(NOT GPU_ARCH)
    message(FATAL_ERROR "GPU_TYPE (${GPU_TYPE}) is not recognized. Available options are [a100, v100].")
  endif()
  
endif()

if(ENABLE_OPENACC)
  add_library(micm_openacc)
  add_library(musica::micm_openacc ALIAS micm_openacc)

  set(GPU_FLAGS -gpu=cc${GPU_ARCH},lineinfo -Minfo=accel)
  message(STATUS "GPU Flags: ${OpenACC_CXX_FLAGS} ${GPU_FLAGS}")

  target_compile_options(micm_openacc PRIVATE ${OpenACC_CXX_FLAGS} ${GPU_FLAGS})

  target_link_libraries(micm_openacc 
    PRIVATE micm
  )

  target_link_libraries(micm_openacc 
    PUBLIC OpenACC::OpenACC_CXX
  )
endif()

if(ENABLE_CUDA)
  message(STATUS "GPU architecture ${GPU_ARCH}")
  set(CMAKE_CUDA_ARCHITECTURES ${GPU_ARCH})

  add_library(micm_cuda)
  add_library(musica::micm_cuda ALIAS micm_cuda)
  set(CMAKE_CUDA_ARCHITECTURES 80)
  target_link_libraries(micm_cuda 
    PRIVATE micm
  )

  target_link_libraries(micm_cuda 
    PUBLIC cudart
  )

  target_compile_definitions(micm_cuda
    PUBLIC
      USE_CUDA
  )
endif()

add_subdirectory(process)