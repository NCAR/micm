include(FetchContent)

################################################################################
# Memory check

if(ENABLE_MEMCHECK)
  find_file(MEMCHECK_SUPPRESS_FILE
    DOC "Suppression file for memory checking"
    NAMES openmpi-valgrind.supp
    PATHS
      /usr/share/openmpi
      /usr/lib64/openmpi/share
      /usr/lib64/openmpi/share/openmpi
      /usr/share)
  if(MEMCHECK_SUPPRESS_FILE)
    set(MEMCHECK_SUPPRESS "--suppressions=${PROJECT_SOURCE_DIR}/test/valgrind.supp --suppressions=${MEMCHECK_SUPPRESS_FILE}")
  else()
    set(MEMCHECK_SUPPRESS "--suppressions=${PROJECT_SOURCE_DIR}/test/valgrind.supp")
  endif()
endif()

################################################################################
# OpenMP

if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  message(STATUS "Compiling with OpenMP support")
endif()

################################################################################
# MPI

if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  message(STATUS "Compiling with MPI support")
endif()

################################################################################
# google test

if(PROJECT_IS_TOP_LEVEL AND ENABLE_LOCAL_TESTS)
  # if google test isn't installed, fetch content will download and build what is needed
  # but, we don't want to run clang tidy on google test, save those variables and reset them later
  foreach (lang IN ITEMS C CXX)
    set("CMAKE_${lang}_CLANG_TIDY_save" "${CMAKE_${lang}_CLANG_TIDY}")
    set("CMAKE_${lang}_CLANG_TIDY" "")
  endforeach ()

  FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.8.1
    # FIND_PACKAGE_ARGS GTest
  )

  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(googletest)

  foreach (lang IN ITEMS C CXX)
    set("CMAKE_${lang}_CLANG_TIDY" "${CMAKE_${lang}_CLANG_TIDY_save}")
  endforeach ()
endif()

################################################################################
# nlohmann::json

if(ENABLE_JSON)
  FetchContent_Declare(json
      GIT_REPOSITORY https://github.com/nlohmann/json.git
      GIT_TAG v3.11.2
  )
  FetchContent_MakeAvailable(json)
endif()

################################################################################
# Docs

if(BUILD_DOCS)
  find_package(Doxygen REQUIRED)
  find_package(Sphinx REQUIRED)
endif()