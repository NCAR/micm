################################################################################
# micm version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/micm/version.hpp @ONLY)

################################################################################
# micm solver
add_library(micm INTERFACE)
add_library(musica::micm ALIAS micm)

if(MICM_COMPILE_WARNING_AS_ERROR)
    # /wd5030: silence "attribute not recognized" for [[gnu::always_inline]]
    # /wd4456: silence false-positive variable shadowing in if/else-if chains (std::get_if pattern)
    # /wd4244,4267: narrowing conversions (double/size_t) widespread in template code
    # /wd4127: conditional expression is constant (pre-C++17 if constexpr patterns)
    # /wd4101: unreferenced local variable (template instantiation artifacts)
    # /wd4459: declaration hides global declaration (common in template code)
    # /wd4100: unreferenced formal parameter (policy/tag types)
  target_compile_options(micm INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:/W4;/WX;/wd5030;/wd4456;/wd4244;/wd4267;/wd4127;/wd4101;/wd4459;/wd4100>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:$<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wpedantic;-Werror;-Wno-unused-parameter>>
    # clang-cl detects as CXX_COMPILER_ID=Clang and gets -Wpedantic, enabling many clang-only warnings
    # that are too aggressive for a C++20 template-heavy codebase targeting GCC/MSVC parity.
    $<$<CXX_COMPILER_ID:Clang>:$<$<COMPILE_LANGUAGE:CXX>:
      -Wno-c++98-compat;-Wno-c++98-compat-pedantic;
      -Wno-sign-conversion;
      -Wno-old-style-cast;
      -Wno-documentation;
      -Wno-global-constructors;
      -Wno-implicit-int-float-conversion;
      -Wno-shadow-uncaptured-local;
      -Wno-unsafe-buffer-usage;
      -Wno-unused-exception-parameter;
      -Wno-newline-eof;
      -Wno-shadow;
      -Wno-float-conversion;
      -Wno-covered-switch-default;
      -Wno-deprecated-copy-with-dtor>>
  )
  # GCC 12 has a false-positive -Wrestrict in <bits/char_traits.h> (fixed in GCC 13)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "13")
    target_compile_options(micm INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wno-restrict>)
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
  if(MICM_COMPILE_WARNING_AS_ERROR)
    target_compile_options(micm_cuda PRIVATE
      $<$<COMPILE_LANGUAGE:CUDA>:-Werror>
      $<$<COMPILE_LANGUAGE:CUDA>:all-warnings>
    )
  endif()
  set_property(TARGET micm_cuda PROPERTY CUDA_STANDARD 20)

  add_subdirectory(process)
  add_subdirectory(solver)
  add_subdirectory(util)
endif()


if (MICM_ENABLE_TESTS)
  target_compile_definitions(micm INTERFACE MICM_ENABLE_TESTS)
endif()

