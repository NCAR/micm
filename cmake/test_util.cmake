################################################################################
# Utility functions for creating tests

if(ENABLE_MEMCHECK)
  find_program(MEMORYCHECK_COMMAND "valgrind")
endif()

################################################################################
# build and add a standard test (one linked to the micm library)

function(create_standard_test)
  set(prefix TEST)
  set(singleValues NAME WORKING_DIRECTORY)
  set(multiValues SOURCES LIBRARIES)

  include(CMakeParseArguments)
  cmake_parse_arguments(${prefix} " " "${singleValues}" "${multiValues}" ${ARGN})

  add_executable(test_${TEST_NAME} ${TEST_SOURCES})

  target_link_libraries(test_${TEST_NAME} PUBLIC musica::micm GTest::gtest_main)

  # link additional libraries
  foreach(library ${TEST_LIBRARIES})
    target_link_libraries(test_${TEST_NAME} PUBLIC ${library})
  endforeach()

  if(ENABLE_JSON)
    target_link_libraries(test_${TEST_NAME} PRIVATE nlohmann_json::nlohmann_json)
  endif()

  if(ENABLE_LLVM)
    target_link_libraries(test_${TEST_NAME} PRIVATE ${llvm_libs})
  endif()

  if(NOT DEFINED TEST_WORKING_DIRECTORY)
    set(TEST_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  endif()

  add_micm_test(${TEST_NAME} test_${TEST_NAME} "" ${TEST_WORKING_DIRECTORY})
endfunction(create_standard_test)

################################################################################
# Add a test

function(add_micm_test test_name test_binary test_args working_dir)
  if(ENABLE_MPI)
    add_test(NAME ${test_name}
      COMMAND mpirun -v -np 2 ${CMAKE_BINARY_DIR}/${test_binary} ${test_args}
             WORKING_DIRECTORY ${working_dir})
  else()
    add_test(NAME ${test_name}
             COMMAND ${test_binary} ${test_args}
             WORKING_DIRECTORY ${working_dir})
  endif()
  set(MEMORYCHECK_COMMAND_OPTIONS "--error-exitcode=1 --trace-children=yes --leak-check=full --gen-suppressions=all ${MEMCHECK_SUPPRESS}")
  set(memcheck "${MEMORYCHECK_COMMAND} ${MEMORYCHECK_COMMAND_OPTIONS}")
  separate_arguments(memcheck)
  if(ENABLE_MPI AND MEMORYCHECK_COMMAND AND ENABLE_MEMCHECK)
    add_test(NAME memcheck_${test_name}
      COMMAND mpirun -v -np 2 ${memcheck} ${CMAKE_BINARY_DIR}/${test_binary} ${test_args}
             WORKING_DIRECTORY ${working_dir})
  elseif(MEMORYCHECK_COMMAND AND ENABLE_MEMCHECK)
    add_test(NAME memcheck_${test_name}
             COMMAND ${memcheck} ${CMAKE_BINARY_DIR}/${test_binary} ${test_args}
             WORKING_DIRECTORY ${working_dir})
  endif()
endfunction(add_micm_test)

################################################################################
