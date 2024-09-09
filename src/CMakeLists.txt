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


if(MICM_ENABLE_CONFIG_READER)
  target_link_libraries(micm INTERFACE yaml-cpp::yaml-cpp)
  target_compile_definitions(micm INTERFACE USE_YAML)
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

if(MICM_ENABLE_CUDA)
  add_library(micm_cuda)
  add_library(musica::micm_cuda ALIAS micm_cuda)
  target_link_libraries(micm_cuda
    PRIVATE micm
  )
  target_link_libraries(micm_cuda 
    PUBLIC CUDA::cudart CUDA::cublas
  )
  set_property(TARGET micm_cuda PROPERTY CUDA_STANDARD 20)
endif()

if (MICM_ENABLE_PROFILE)
  target_compile_definitions(micm INTERFACE MICM_PROFILE)
endif()

add_subdirectory(process)
add_subdirectory(solver)
add_subdirectory(util)
