# CMake generated Testfile for 
# Source directory: /Users/mthind/Desktop/code/micm/examples
# Build directory: /Users/mthind/Desktop/code/micm/build/examples
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[carbon_bond_5_example]=] "micm" "/Users/mthind/Desktop/code/micm/build/examples/configs/carbon_bond_5" "/Users/mthind/Desktop/code/micm/build/examples/configs/carbon_bond_5/initial_conditions.csv")
set_tests_properties([=[carbon_bond_5_example]=] PROPERTIES  WORKING_DIRECTORY "/Users/mthind/Desktop/code/micm/build" _BACKTRACE_TRIPLES "/Users/mthind/Desktop/code/micm/examples/CMakeLists.txt;12;add_test;/Users/mthind/Desktop/code/micm/examples/CMakeLists.txt;0;")
add_test([=[chapman_example]=] "micm" "/Users/mthind/Desktop/code/micm/build/examples/configs/chapman" "/Users/mthind/Desktop/code/micm/build/examples/configs/chapman/initial_conditions.csv")
set_tests_properties([=[chapman_example]=] PROPERTIES  WORKING_DIRECTORY "/Users/mthind/Desktop/code/micm/build" _BACKTRACE_TRIPLES "/Users/mthind/Desktop/code/micm/examples/CMakeLists.txt;15;add_test;/Users/mthind/Desktop/code/micm/examples/CMakeLists.txt;0;")
add_test([=[robertson_example]=] "micm" "/Users/mthind/Desktop/code/micm/build/examples/configs/robertson" "/Users/mthind/Desktop/code/micm/build/examples/configs/robertson/initial_conditions.csv")
set_tests_properties([=[robertson_example]=] PROPERTIES  WORKING_DIRECTORY "/Users/mthind/Desktop/code/micm/build" _BACKTRACE_TRIPLES "/Users/mthind/Desktop/code/micm/examples/CMakeLists.txt;18;add_test;/Users/mthind/Desktop/code/micm/examples/CMakeLists.txt;0;")
add_test([=[ts1_example]=] "micm" "/Users/mthind/Desktop/code/micm/build/examples/configs/TS1" "/Users/mthind/Desktop/code/micm/build/examples/configs/TS1/initial_conditions.csv")
set_tests_properties([=[ts1_example]=] PROPERTIES  WORKING_DIRECTORY "/Users/mthind/Desktop/code/micm/build" _BACKTRACE_TRIPLES "/Users/mthind/Desktop/code/micm/examples/CMakeLists.txt;21;add_test;/Users/mthind/Desktop/code/micm/examples/CMakeLists.txt;0;")
