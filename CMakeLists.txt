################################################################################
# Preamble
cmake_minimum_required(VERSION 3.21)

project(
  micm
  VERSION 3.0.0
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

include(GNUInstallDirs)
set(INSTALL_PREFIX "micm-${PROJECT_VERSION}" )
set(INSTALL_MOD_DIR "${INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")

option(ENABLE_CLANG_TIDY "Enable clang-tiday to format source code" OFF)
option(ENABLE_MPI "Enable MPI parallel support" OFF)
option(ENABLE_OPENMP "Enable OpenMP support" OFF)
option(ENABLE_COVERAGE "Enable code coverage output" OFF)
option(ENABLE_MEMCHECK "Enable memory checking in tests" OFF)
option(ENABLE_JSON "Enable json configureation file reading" ON)
option(ENABLE_REGRESSION_TESTS "Enable regression tests against the old pre-processed version of micm" ON)
option(BUILD_DOCS "Build the documentation" OFF)

################################################################################
# Dependencies

include(cmake/dependencies.cmake)

################################################################################
# micm targets and documentation

add_subdirectory(src)

if(BUILD_DOCS)
  add_subdirectory(docs)
endif()

################################################################################
# Tests

if(PROJECT_IS_TOP_LEVEL)
  # Test code coverage
  if(ENABLE_COVERAGE)
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
