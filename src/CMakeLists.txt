################################################################################
# micm version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/micm/version.hpp @ONLY)

################################################################################
# micm solver
add_library(micm INTERFACE)
add_library(musica::micm ALIAS micm)

target_compile_features(micm INTERFACE cxx_std_20)
set_property(TARGET micm PROPERTY INTERFACE_CXX_STANDARD 20)
set_property(TARGET micm PROPERTY INTERFACE_CXX_STANDARD_REQUIRED YES)

target_compile_definitions(micm INTERFACE MICM_DEFAULT_VECTOR_SIZE=${MICM_DEFAULT_VECTOR_SIZE})

# Propagate the precision both ways.  include/micm/util/types.hpp switches on MICM_USE_SINGLE so
# that double precision stays the default for anyone including the headers without this target;
# MICM_USE_DOUBLE is still defined when the option is ON so consumers can test for either.
if(MICM_USE_DOUBLE)
  target_compile_definitions(micm INTERFACE MICM_USE_DOUBLE)
else()
  target_compile_definitions(micm INTERFACE MICM_USE_SINGLE)
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
  
  # --expt-relaxed-constexpr is nvcc-only; guard it so non-nvcc CUDA compilers
  # (e.g. clang, used by the clang-tidy CI to generate parseable compile commands)
  # don't choke on it.
  target_compile_options(micm_cuda PRIVATE $<$<CUDA_COMPILER_ID:NVIDIA>:--expt-relaxed-constexpr>)
  set_property(TARGET micm_cuda PROPERTY CUDA_STANDARD 20)

  add_subdirectory(process)
  add_subdirectory(solver)
  add_subdirectory(util)
endif()

if(MICM_ENABLE_KOKKOS)
  # MICM_USE_KOKKOS switches MICM_LAMBDA, MICM_DEVICE_FUNCTION and MICM_CONSTEXPR for
  # every micm header, so it has to be one state for the whole build.  Defining it only
  # on micm_kokkos while micm also linked Kokkos let two translation units in the same
  # binary compile the same header-only functions with different macro expansions --
  # an ODR violation.  Carry both the link and the definition on micm itself.
  target_link_libraries(micm INTERFACE Kokkos::kokkos)
  target_compile_definitions(micm INTERFACE MICM_USE_KOKKOS)

  if(NOT MICM_KOKKOS_DEFAULT_VECTOR_SIZE STREQUAL "")
    target_compile_definitions(micm INTERFACE MICM_KOKKOS_DEFAULT_VECTOR_SIZE=${MICM_KOKKOS_DEFAULT_VECTOR_SIZE})
  endif()

  # Retained so existing `musica::micm_kokkos` consumers keep working; micm already
  # provides Kokkos and the macro, so this adds nothing beyond the name.
  add_library(micm_kokkos INTERFACE)
  add_library(musica::micm_kokkos ALIAS micm_kokkos)
  target_link_libraries(micm_kokkos INTERFACE micm)
endif()

if (MICM_ENABLE_TESTS)
  target_compile_definitions(micm INTERFACE MICM_ENABLE_TESTS)
endif()

