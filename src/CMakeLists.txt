################################################################################
# micm version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/micm/version.hpp @ONLY)

################################################################################
# micm solver
add_library(micm INTERFACE)
add_library(musica::micm ALIAS micm)

target_compile_features(micm INTERFACE cxx_std_20)

target_compile_definitions(micm INTERFACE MICM_DEFAULT_VECTOR_SIZE=${MICM_DEFAULT_VECTOR_MATRIX_SIZE})

target_include_directories(micm
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Check for Intel compiler
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel" OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  # Set compiler flags for Intel C++
  if (WIN32)
    # For Windows
    target_compile_options(micm INTERFACE "/fp:precise")
  else ()
    # For other operating systems (e.g., Linux, macOS)
    target_compile_options(micm INTERFACE "-ffp-model=precise")
  endif ()
endif ()


if(MICM_ENABLE_JSON)
  target_link_libraries(micm INTERFACE nlohmann_json::nlohmann_json)
  target_compile_definitions(micm INTERFACE USE_JSON)
endif()

if(MICM_ENABLE_OPENMP)
  target_link_libraries(micm INTERFACE OpenMP::OpenMP_CXX)
endif()

if(MICM_ENABLE_MPI)
  target_link_libraries(micm INTERFACE MPI::MPI_CXX)
endif()

if(MICM_ENABLE_LLVM)
  target_link_libraries(micm INTERFACE ${llvm_libs})
endif()

if(MICM_ENABLE_OPENACC OR MICM_ENABLE_CUDA)
  if(NOT MICM_GPU_TYPE)
    message(FATAL_ERROR "MICM_GPU_TYPE is not set or is empty. Please provide a GPU type.")
  endif()

  set(MICM_GPU_ARCH "")
  if (MICM_GPU_TYPE STREQUAL a100)
    set(MICM_GPU_ARCH "80")
  endif()
  if (MICM_GPU_TYPE STREQUAL v100)
    set(MICM_GPU_ARCH "70")
  endif()

  if(NOT MICM_GPU_ARCH)
    message(FATAL_ERROR "MICM_GPU_TYPE (${MICM_GPU_TYPE}) is not recognized. Available options are [a100, v100].")
  endif()
  
endif()

if(MICM_ENABLE_OPENACC)
  add_library(micm_openacc)
  add_library(musica::micm_openacc ALIAS micm_openacc)

  set(GPU_FLAGS -gpu=cc${MICM_GPU_ARCH},lineinfo -Minfo=accel)
  message(STATUS "GPU Flags: ${OpenACC_CXX_FLAGS} ${GPU_FLAGS}")

  target_compile_options(micm_openacc PRIVATE ${OpenACC_CXX_FLAGS} ${GPU_FLAGS})

  target_link_libraries(micm_openacc 
    PRIVATE micm
  )

  target_link_libraries(micm_openacc 
    PUBLIC OpenACC::OpenACC_CXX
  )
endif()

if(MICM_ENABLE_CUDA)
  message(STATUS "GPU architecture ${MICM_GPU_ARCH}")
  set(CMAKE_CUDA_ARCHITECTURES ${MICM_GPU_ARCH})

  add_library(micm_cuda)
  add_library(musica::micm_cuda ALIAS micm_cuda)
  set(CMAKE_CUDA_ARCHITECTURES 80)
  target_link_libraries(micm_cuda 
    PRIVATE micm
  )

  target_link_libraries(micm_cuda 
    PUBLIC cudart cublas culibos
  )

endif()

add_subdirectory(process)
add_subdirectory(solver)
add_subdirectory(util)
