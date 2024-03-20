################################################################################
# Preamble
cmake_minimum_required(VERSION 3.21)

project(
  micm
  VERSION 3.4.0
  LANGUAGES CXX
)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

message(STATUS "CMake build configuration for micm(${CMAKE_BUILD_TYPE}) ${PROJECT_VERSION}")

################################################################################
# Projet wide setup options

include(cmake/StaticAnalyzers.cmake)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")

# Set up include and lib directories
set(MICM_LIB_DIR "${PROJECT_BINARY_DIR}/lib")

option(MICM_ENABLE_CLANG_TIDY "Enable clang-tiday to format source code" OFF)
option(MICM_ENABLE_MPI "Enable MPI parallel support" OFF)
option(MICM_ENABLE_OPENMP "Enable OpenMP support" OFF)
option(MICM_ENABLE_COVERAGE "Enable code coverage output" OFF)
option(MICM_ENABLE_MEMCHECK "Enable memory checking in tests" OFF)
option(MICM_ENABLE_JSON "Enable json configureation file reading" ON)
option(MICM_BUILD_DOCS "Build the documentation" OFF)
option(MICM_ENABLE_CUDA "Build with Cuda support" OFF)
option(MICM_ENABLE_OPENACC "Build with OpenACC Support" OFF)
option(MICM_ENABLE_LLVM "Build with LLVM support for JIT-compiling" OFF)
option(MICM_ENABLE_TESTS "Build the tests" ON)
option(MICM_ENABLE_EXAMPLES "Build the examples" ON)
set(MICM_DEFAULT_VECTOR_MATRIX_SIZE "4" CACHE STRING "Default size for vectorizable matrix types")

include(CMakeDependentOption)
# Option to collect custom OpenACC flags
cmake_dependent_option(MICM_ENABLE_GPU_TYPE "Enable custom CUDA flags" OFF
  "MICM_ENABLE_OPENACC OR MICM_ENABLE_CUDA" ON)

if(MICM_ENABLE_GPU_TYPE)
  set(MICM_GPU_TYPE "" CACHE STRING "The GPU type being targeted")
endif()

# on ubuntu with clang, an incorrect version of the c++ standard library was being linked
if (${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux" AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # If the compiler is Clang, use libc++
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

################################################################################
# Dependencies

include(cmake/dependencies.cmake)

################################################################################
# micm targets and documentation

add_subdirectory(src)

if(MICM_BUILD_DOCS)
  add_subdirectory(docs)
endif()

################################################################################
# Tests

if(PROJECT_IS_TOP_LEVEL AND MICM_ENABLE_TESTS)
  # Test code coverage
  if(MICM_ENABLE_COVERAGE)
    include(CodeCoverage)
    append_coverage_compiler_flags()
    setup_target_for_coverage_lcov(
        NAME coverage
        EXECUTABLE "ctest"
        EXCLUDE "${PROJECT_SOURCE_DIR}/test/*"
        BASE_DIRECTORY "${PROJECT_SOURCE_DIR}/src")

  endif()

  enable_testing()
  add_subdirectory(test)
endif()

################################################################################
# Packaging

# only include packaging if we are the top level project being built
if(PROJECT_IS_TOP_LEVEL)
  add_subdirectory(packaging)
endif()

################################################################################

# on ubuntu with clang, an incorrect version of the c++ standard library was being linked
if (${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux" AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # If the compiler is Clang, use libc++
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

################################################################################
# Command-line driver examples

if(PROJECT_IS_TOP_LEVEL AND MICM_ENABLE_EXAMPLES AND MICM_ENABLE_JSON)
  add_subdirectory(examples)
endif()

################################################################################
# Copy configuration data

if(PROJECT_IS_TOP_LEVEL AND (MICM_ENABLE_TESTS OR MICM_ENABLE_EXAMPLES))
  add_custom_target(copy_example_configs ALL ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/examples/configs ${CMAKE_BINARY_DIR}/examples/configs)
endif()
