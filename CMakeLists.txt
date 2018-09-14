# Set default project to unknown
if(NOT PROJECT)
  message(STATUS "Setting CCPP project to 'unknown' as none was specified.")
  set(PROJECT "Unknown")
endif (NOT PROJECT)

#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.11)

# Use rpaths on MacOSX
set(CMAKE_MACOSX_RPATH 1)

if(POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
    project(micmchem VERSION 1.0.0)
else(POLICY CMP0048)
    project(micmchem)
    set(PROJECT_VERSION 1.0.0)
    set(PROJECT_VERSION_MAJOR 1)
    set(PROJECT_VERSION_MINOR 0)
    set(PROJECT_VERSION_PATCH 0)
endif(POLICY CMP0048)

if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif(POLICY CMP0042)

#------------------------------------------------------------------------------
set(PACKAGE "micm-chemistry")
set(AUTHORS  "Cheryl Craig" "Francis Vitt" "Andrew Conley" "Stacy Walters")

#------------------------------------------------------------------------------
# Enable Fortran
enable_language(Fortran)

# DH this should be a PROJECT == CCPP-SCM test
if (NOT PROJECT STREQUAL "CCPP-FV3")
  #------------------------------------------------------------------------------
  # CMake Modules
  # Set the CMake module path
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../ccpp-framework/cmake")
endif (NOT PROJECT STREQUAL "CCPP-FV3")

#------------------------------------------------------------------------------
# Set OpenMP flags for C/C++/Fortran
if (OPENMP)
  include(detect_openmp)
  detect_openmp()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  message(STATUS "Enable OpenMP support for C/C++/Fortran compiler")
else(OPENMP)
  message (STATUS "Disable OpenMP support for C/C++/Fortran compiler")
endif()

#------------------------------------------------------------------------------
# The Fortran compiler/linker flag inserted by cmake to create shared libraries
# with the Intel compiler is deprecated (-i_dynamic), correct here.
# CMAKE_Fortran_COMPILER_ID = {"Intel", "PGI", "GNU", "Clang", "MSVC", ...}
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    string(REPLACE "-i_dynamic" "-shared-intel"
           CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS
           "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS}")
    string(REPLACE "-i_dynamic" "-shared-intel"
           CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS
           "${CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS}")
endif()

#------------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "Coverage")
endif()

#------------------------------------------------------------------------------
# By default we want a shared library
option(BUILD_SHARED_LIBS "Build a shared library" ON)

#------------------------------------------------------------------------------
# Add the CCPP include/module directory
set(CCPP_INCLUDE_DIRS "" CACHE FILEPATH "Path to ccpp includes")
set_property(DIRECTORY PROPERTY INCLUDE_DIRECTORIES ${CCPP_INCLUDE_DIRS})

#------------------------------------------------------------------------------
# Add the CCPP library
set(CCPP_LIB_DIRS "" CACHE FILEPATH "Path to ccpp library")
link_directories(${CCPP_LIB_DIRS})
list(APPEND LIBS "ccpp")

# Add the required preprocessor flags so that cmake can sort out the dependencies
ADD_DEFINITIONS(-DNEMS_GSM)

#------------------------------------------------------------------------------
# Set the sources: physics schemes
if (PROJECT STREQUAL "CCPP-FV3")
  include(./CCPP_SCHEMES.cmake)
else (PROJECT STREQUAL "CCPP-FV3")
  set(SCHEMES
      ./src/machine.F
      ./src/chem_solve.F90
      ./src/kinetics.F90
      ./src/k_rateConst_terminator.F90
      ./src/ODE_solver.F90
      ./src/kinetics_module.F90
      ./src/half_solver.F90
      ./src/rosenbrock_solver.F90
      ./src/mozart_solver.F90
      ./src/solver_var_defs.f90
  )
endif (PROJECT STREQUAL "CCPP-FV3")

# Set the sources: physics scheme caps
set(CAPS
      ./caps/chem_solve_cap.F90
      ./caps/kinetics_cap.F90
      ./caps/k_rateConst_terminator_cap.F90
  )

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
else (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${CMAKE_Fortran_COMPILER_ID})
  message (FATAL_ERROR "This program has only been compiled with gfortran, pgf90 and ifort. If another compiler is needed, the appropriate flags must be added in ${GFS_PHYS_SRC}/CMakeLists.txt")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

#apply general fortran tags to all fortran source files
if(${CMAKE_VERSION} LESS 3.3)
  string (REPLACE ";" " " f_flags_str "${f_flags}")
  SET_PROPERTY(SOURCE ${SCHEMES} ${CAPS} APPEND_STRING PROPERTY COMPILE_FLAGS " ${f_flags_str}")
else(${CMAKE_VERSION} LESS 3.3)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${f_flags}>")
endif (${CMAKE_VERSION} LESS 3.3)

#if (PROJECT STREQUAL "CCPP-FV3")
#  link_directories(${NCEPLIBS_DIR}/lib)
## DH this should be a PROJECT == CCPP-SCM and else undefined test (similar in other places)
#else (PROJECT STREQUAL "CCPP-FV3")
#  SET(W3LIB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../external/w3nco/v2.0.6/src)
#  SET(BACIOLIB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../external/bacio/v2.0.1/src)
#  SET(SPLIB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../external/sp/v2.0.2/src)
#
#  #add "sibling" directories (must specify the build directory too)
#  ADD_SUBDIRECTORY(${W3LIB_SRC} ${CMAKE_BINARY_DIR}/w3nco)
#  ADD_SUBDIRECTORY(${BACIOLIB_SRC} ${CMAKE_BINARY_DIR}/bacio)
#  ADD_SUBDIRECTORY(${SPLIB_SRC} ${CMAKE_BINARY_DIR}/sp)
#
#  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/w3nco)
#  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/sp)
#  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/bacio)
#
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/ccpp-framework/src)
#endif (PROJECT STREQUAL "CCPP-FV3")

#------------------------------------------------------------------------------

add_library(micmchem ${SCHEMES} ${CAPS})
if (PROJECT STREQUAL "CCPP-FV3")
  target_link_libraries(micmchem LINK_PUBLIC ${LIBS} sp_v2.0.2_d bacio_4 w3nco_d)
else (PROJECT STREQUAL "CCPP-FV3")
  target_link_libraries(micmchem LINK_PUBLIC ${LIBS} )
endif (PROJECT STREQUAL "CCPP-FV3")
set_target_properties(micmchem PROPERTIES VERSION ${PROJECT_VERSION}
                                     SOVERSION ${PROJECT_VERSION_MAJOR}
                                     COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                                     LINK_FLAGS "${CMAKE_Fortran_FLAGS}")
# DH* hack for PGI compiler: rename objects in scheme cap object files for ISO_C compliancy
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
set(CAPOBJS)
foreach(cap ${CAPS})
    string(REPLACE "_cap.F90" "_cap.F90.o" capobj "./${CMAKE_FILES_DIRECTORY}/micmchem.dir/${cap}")
    list(APPEND CAPOBJS ${capobj})
endforeach(cap)

add_custom_command(TARGET micmchem
                   PRE_LINK
                   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/pgifix.py --cmake ${CAPOBJS}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   COMMENT "Running pgifix_wrapper.py over all scheme caps")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
# *DH end hack for PGI compiler

if (PROJECT STREQUAL "CCPP-FV3")
  # Define where to install the library
  install(TARGETS micmchem
          EXPORT micmchem-targets
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib
          RUNTIME DESTINATION lib
  )
  # Export our configuration
  install(EXPORT micmchem-targets
          FILE micmchem-config.cmake
          DESTINATION lib/cmake
  )
endif (PROJECT STREQUAL "CCPP-FV3")
