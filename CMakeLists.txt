################################################################################
# Preamble
cmake_minimum_required(VERSION 3.21)

project(
  micm
  VERSION 3.0.0
  LANGUAGES CXX
)

message(STATUS "CMake build configuration for micm(${CMAKE_BUILD_TYPE}) ${PROJECT_VERSION}")

################################################################################
# Projet wide setup options

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")

# Set up include and lib directories
set(MICM_LIB_DIR "${PROJECT_BINARY_DIR}/lib")

include(GNUInstallDirs)
set(
  INSTALL_PREFIX 
  "${CMAKE_INSTALL_PREFIX}/ncar/${PROJECT_NAME}-${PROJECT_VERSION}"
)

option(ENABLE_MPI "Enable MPI parallel support" OFF)
option(ENABLE_OPENMP "Enable OpenMP support" OFF)
option(ENABLE_COVERAGE "Enable code coverage output" OFF)
option(ENABLE_MEMCHECK "Enable memory checking in tests" OFF)

################################################################################
# Dependencies

include(cmake/dependencies.cmake)

################################################################################
# micm targets

add_subdirectory(src)

################################################################################
# Tests

if(PROJECT_IS_TOP_LEVEL)
  # Test code coverage
  if(ENABLE_COVERAGE)
    include(CodeCoverage)
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
