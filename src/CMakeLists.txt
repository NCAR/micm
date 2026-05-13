################################################################################
# micm version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/micm/version.hpp @ONLY)

################################################################################
# micm solver
add_library(micm INTERFACE)
add_library(musica::micm ALIAS micm)

if(MICM_COMPILE_WARNING_AS_ERROR)
  if(MSVC)
    # MSVC is true for both MSVC and clang-cl (MSVC frontend)
    # /wd5030: silence "attribute not recognized" for [[gnu::always_inline]]
    # /wd4456: silence false-positive variable shadowing in if/else-if chains (std::get_if pattern)
    target_compile_options(micm INTERFACE /W4 /WX /wd5030 /wd4456)
  else()
    target_compile_options(micm INTERFACE -Wall -Werror)
    # GCC 12 has a false-positive -Wrestrict in <bits/char_traits.h> (fixed in GCC 13)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "13")
      target_compile_options(micm INTERFACE -Wno-restrict)
    endif()
  endif()
endif()

target_compile_features(micm INTERFACE cxx_std_20)
set_property(TARGET micm PROPERTY INTERFACE_CXX_STANDARD 20)
set_property(TARGET micm PROPERTY INTERFACE_CXX_STANDARD_REQUIRED YES)

target_compile_definitions(micm INTERFACE MICM_DEFAULT_VECTOR_SIZE=${MICM_DEFAULT_VECTOR_SIZE})

# MSVC's <cmath> requires _USE_MATH_DEFINES to expose M_PI etc.
# Define it as a compile definition rather than in source to avoid the reserved-identifier warning.
if(MSVC)
  target_compile_definitions(micm INTERFACE _USE_MATH_DEFINES)
endif()

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

if(MICM_ENABLE_OPENMP)
  target_link_libraries(micm INTERFACE OpenMP::OpenMP_CXX)
endif()

if(MICM_ENABLE_MPI)
  target_link_libraries(micm INTERFACE MPI::MPI_CXX)
endif()

if(NOT ${MICM_GPU_TYPE} STREQUAL "None")
  if(MICM_BUILD_SHARED_LIBS)
    add_library(micm_cuda SHARED)
  else()
    add_library(micm_cuda STATIC)
  endif()

  add_library(musica::micm_cuda ALIAS micm_cuda)
  target_link_libraries(micm_cuda
    PRIVATE micm
  )
  
  target_link_libraries(micm_cuda 
    PUBLIC CUDA::cudart CUDA::cublas
  )
  
  target_compile_options(micm_cuda PRIVATE --expt-relaxed-constexpr)
  set_property(TARGET micm_cuda PROPERTY CUDA_STANDARD 20)

  add_subdirectory(process)
  add_subdirectory(solver)
  add_subdirectory(util)
endif()


if (MICM_ENABLE_TESTS)
  target_compile_definitions(micm INTERFACE MICM_ENABLE_TESTS)
endif()

