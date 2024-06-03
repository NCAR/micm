include(FetchContent)

################################################################################
# Memory check

if(MICM_ENABLE_MEMCHECK)
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

if(MICM_ENABLE_OPENMP)
  if(APPLE)
    # Apple clang by default doesn't include support for openmp
    # but if omp was installed with `brew install libomp`, support can be configured
    if(CMAKE_C_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      # Set the C flags
      set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
      set(OpenMP_C_LIB_NAMES "omp")
      set(OpenMP_omp_LIBRARY omp)

      # Set the CXX flags
      set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
      set(OpenMP_CXX_LIB_NAMES "omp")
      set(OpenMP_omp_LIBRARY omp)

      # Assume that libomp is instaleld on mac with brew when using apple clang
      # Get the path to libomp from Homebrew
      execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE LIBOMP_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )

      # Set the full path to the libomp library
      set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")

      # Set the include directory
      set(LIBOMP_INCLUDE_DIR "${LIBOMP_PREFIX}/include")

      include_directories(
        ${LIBOMP_INCLUDE_DIR}
      )
    endif()

  endif()

  find_package(OpenMP REQUIRED)
  message(STATUS "Compiling with OpenMP support")
endif()

################################################################################
# MPI

if(MICM_ENABLE_MPI)
  find_package(MPI REQUIRED)
  message(STATUS "Compiling with MPI support")
endif()

################################################################################
# google test

if(MICM_ENABLE_TESTS)
  FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG be03d00f5f0cc3a997d1a368bee8a1fe93651f48
  )

  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(googletest)

  # don't run clang-tidy on google test
  set_target_properties(gtest PROPERTIES CXX_CLANG_TIDY "")
  set_target_properties(gtest_main PROPERTIES CXX_CLANG_TIDY "")
  # set_target_properties(gmock PROPERTIES CXX_CLANG_TIDY "")
  # set_target_properties(gmock_main PROPERTIES CXX_CLANG_TIDY "")
endif()

################################################################################
# nlohmann::json

if(MICM_ENABLE_JSON)
  FetchContent_Declare(json
      GIT_REPOSITORY https://github.com/nlohmann/json.git
      GIT_TAG v3.11.2
  )
  FetchContent_MakeAvailable(json)
endif()

################################################################################
# Docs

if(MICM_BUILD_DOCS)
  find_package(Doxygen REQUIRED)
  find_package(Sphinx REQUIRED)
endif()

################################################################################
# GPU Support

if(MICM_ENABLE_CUDA)
  if(NOT MICM_GPU_TYPE)
    message(FATAL_ERROR "MICM_GPU_TYPE is not set or is empty. Please provide a GPU type.")
  endif()

  set(MICM_GPU_ARCH "")
  if (MICM_GPU_TYPE STREQUAL a100)
  set(MICM_GPU_ARCH "80")
  endif()
  if (MICM_GPU_TYPE STREQUAL v100)
    set(MICM_GPU_ARCH "70")
  endif()
  if(NOT MICM_GPU_ARCH)
    message(FATAL_ERROR "MICM_GPU_TYPE (${MICM_GPU_TYPE}) is not recognized. Available options are [a100, v100].")
  endif()

  message(STATUS "GPU architecture ${MICM_GPU_ARCH}")

  include(CheckLanguage)
  check_language(CUDA)

  message(STATUS "CMAKE_CUDA_COMPILER = ${CMAKE_CUDA_COMPILER}")

  if(NOT CMAKE_CUDA_COMPILER)
    message(FATAL_ERROR "Unable to find compatiable compiler for CUDA.")
  endif()

  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
endif()

################################################################################
# LLVM Support
#
# TODO: Try to use fetch content for LLVM libraries

if(MICM_ENABLE_LLVM)
  find_package(LLVM REQUIRED CONFIG)
  if(LLVM_FOUND)
    message(STATUS "Found LLVM ${LLVM_PACKAGE_VERSION}")
    message(STATUS "Using LLVMConfig.cmake in: ${LLVM_DIR}")

    include_directories(${LLVM_INCLUDE_DIRS})
    separate_arguments(LLVM_DEFINITIONS_LIST NATIVE_COMMAND ${LLVM_DEFINITIONS})
    add_definitions(${LLVM_DEFINITIONS_LIST})

    llvm_map_components_to_libnames(llvm_libs support core orcjit native irreader)
  else()
    set(LLVM_CMD "llvm-config --cxxflags --ldflags --system-libs --libs support core orcjit native irreader | tr '\\n' ' '")
    execute_process(COMMAND bash "-c" ${LLVM_CMD}
                    OUTPUT_VARIABLE llvm_libs)
    separate_arguments(llvm_libs)
  endif()
endif()
