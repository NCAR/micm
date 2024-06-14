# disable clang tidy for test files
set(CMAKE_CXX_CLANG_TIDY "")

if(MICM_ENABLE_EXAMPLES)
  add_executable(micmDriver example.cpp util_example.cpp)
  target_link_libraries(micmDriver PUBLIC musica::micm)
  set_target_properties(micmDriver PROPERTIES OUTPUT_NAME "micm")

  ################################################################################
  # Run each example configuration as a test

  add_test(NAME carbon_bond_5_example
          COMMAND micm ${CMAKE_BINARY_DIR}/examples/configs/carbon_bond_5 ${CMAKE_BINARY_DIR}/examples/configs/carbon_bond_5/initial_conditions.csv
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  add_test(NAME chapman_example
          COMMAND micm ${CMAKE_BINARY_DIR}/examples/configs/chapman ${CMAKE_BINARY_DIR}/examples/configs/chapman/initial_conditions.csv
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  add_test(NAME robertson_example
          COMMAND micm ${CMAKE_BINARY_DIR}/examples/configs/robertson ${CMAKE_BINARY_DIR}/examples/configs/robertson/initial_conditions.csv
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  add_test(NAME ts1_example
          COMMAND micm ${CMAKE_BINARY_DIR}/examples/configs/TS1 ${CMAKE_BINARY_DIR}/examples/configs/TS1/initial_conditions.csv
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

if(MICM_ENABLE_PROFILE)
  add_executable(micmprofiler profile_example.cpp)
  target_link_libraries(micmprofiler PUBLIC musica::micm)
  set_target_properties(micmprofiler PROPERTIES OUTPUT_NAME "profile")

  # Run TS1 configuration as a test
  add_test(NAME profile_ts1
          COMMAND profile ${CMAKE_BINARY_DIR}/examples/configs/TS1 ${CMAKE_BINARY_DIR}/examples/configs/TS1/initial_conditions.csv
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()