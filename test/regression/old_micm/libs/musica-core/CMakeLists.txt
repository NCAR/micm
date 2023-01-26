cmake_minimum_required(VERSION 3.12)
set(CMAKE_USER_MAKE_RULES_OVERRIDE "SetDefaults.cmake")
project(MUSICA)
enable_language(Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-modules)

################################################################################
# options

option(ENABLE_COVERAGE "Enable code coverage output" ON)

################################################################################
# Test code coverage

if(ENABLE_COVERAGE)
  include(CodeCoverage)
  if(CMAKE_COMPILER_IS_GNUCXX)
    setup_target_for_coverage_lcov(
        NAME coverage
        EXECUTABLE "ctest"
        EXCLUDE "${PROJECT_SOURCE_DIR}/test/*"
        BASE_DIRECTORY "${PROJECT_SOURCE_DIR}/src")
  endif()
endif()

################################################################################
# NetCDF library

find_path(NETCDF_INCLUDE_DIR netcdf.mod NETCDF.mod
  DOC "NetCDF include directory (must contain netcdf.mod)"
  PATHS
    $ENV{NETCDF_HOME}/include
    /usr/lib/gfortran/modules
    /usr/lib64/gfortran/modules
    /opt/local/include)
find_library(NETCDF_C_LIB netcdf
  DOC "NetCDF C library"
  PATHS
    $ENV{NETCDF_HOME}/lib
    $ENV{NETCDF_HOME}/lib64
    opt/local/lib)
find_library(NETCDF_FORTRAN_LIB netcdff
  DOC "NetCDF Fortran library"
  PATHS
    $ENV{NETCDF_HOME}/lib
    $ENV{NETCDF_HOME}/lib64
    /opt/local/lib)
set(NETCDF_LIBS ${NETCDF_C_LIB})
if(NETCDF_FORTRAN_LIB)
  set(NETCDF_LIBS ${NETCDF_LIBS} ${NETCDF_FORTRAN_LIB})
endif()
include_directories(${NETCDF_INCLUDE_DIR})

################################################################################
# json-fortran library

find_path(JSON_INCLUDE_DIR json_module.mod
  DOC "json-fortran include directory (must include json_*.mod files)"
  PATHS
    $ENV{JSON_FORTRAN_HOME}/lib
    /opt/local/lib
    /usr/local/lib
    /usr/local/lib64)
find_library(JSON_LIB jsonfortran
  DOC "json-fortran library"
  PATHS
    $ENV{JSON_FORTRAN_HOME}/lib
    /opt/local/lib
    /usr/local/lib
    /usr/local/lib64)
include_directories(${JSON_INCLUDE_DIR})

################################################################################
# Add flags for various compilers

if(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  add_definitions(-DMUSICA_USING_INTEL)
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  add_definitions(-DMUSICA_USING_GNU)
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
  add_definitions(-DMUSICA_USING_PGI)
endif()

################################################################################
# MUSICA library

add_subdirectory(src)
include_directories(${CMAKE_BINARY_DIR}/src)

################################################################################
# MUSICA tests

enable_testing()
add_subdirectory(test)
