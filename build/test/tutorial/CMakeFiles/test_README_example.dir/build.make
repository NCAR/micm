# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.30.2/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.30.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/mthind/Desktop/code/micm

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/mthind/Desktop/code/micm/build

# Include any dependencies generated for this target.
include test/tutorial/CMakeFiles/test_README_example.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/tutorial/CMakeFiles/test_README_example.dir/compiler_depend.make

# Include the progress variables for this target.
include test/tutorial/CMakeFiles/test_README_example.dir/progress.make

# Include the compile flags for this target's objects.
include test/tutorial/CMakeFiles/test_README_example.dir/flags.make

test/tutorial/CMakeFiles/test_README_example.dir/test_README_example.cpp.o: test/tutorial/CMakeFiles/test_README_example.dir/flags.make
test/tutorial/CMakeFiles/test_README_example.dir/test_README_example.cpp.o: /Users/mthind/Desktop/code/micm/test/tutorial/test_README_example.cpp
test/tutorial/CMakeFiles/test_README_example.dir/test_README_example.cpp.o: test/tutorial/CMakeFiles/test_README_example.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/mthind/Desktop/code/micm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/tutorial/CMakeFiles/test_README_example.dir/test_README_example.cpp.o"
	cd /Users/mthind/Desktop/code/micm/build/test/tutorial && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/tutorial/CMakeFiles/test_README_example.dir/test_README_example.cpp.o -MF CMakeFiles/test_README_example.dir/test_README_example.cpp.o.d -o CMakeFiles/test_README_example.dir/test_README_example.cpp.o -c /Users/mthind/Desktop/code/micm/test/tutorial/test_README_example.cpp

test/tutorial/CMakeFiles/test_README_example.dir/test_README_example.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test_README_example.dir/test_README_example.cpp.i"
	cd /Users/mthind/Desktop/code/micm/build/test/tutorial && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mthind/Desktop/code/micm/test/tutorial/test_README_example.cpp > CMakeFiles/test_README_example.dir/test_README_example.cpp.i

test/tutorial/CMakeFiles/test_README_example.dir/test_README_example.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test_README_example.dir/test_README_example.cpp.s"
	cd /Users/mthind/Desktop/code/micm/build/test/tutorial && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mthind/Desktop/code/micm/test/tutorial/test_README_example.cpp -o CMakeFiles/test_README_example.dir/test_README_example.cpp.s

# Object files for target test_README_example
test_README_example_OBJECTS = \
"CMakeFiles/test_README_example.dir/test_README_example.cpp.o"

# External object files for target test_README_example
test_README_example_EXTERNAL_OBJECTS =

test_README_example: test/tutorial/CMakeFiles/test_README_example.dir/test_README_example.cpp.o
test_README_example: test/tutorial/CMakeFiles/test_README_example.dir/build.make
test_README_example: lib/libgtest_main.a
test_README_example: _deps/yaml-cpp-build/libyaml-cpp.a
test_README_example: lib/libgtest.a
test_README_example: test/tutorial/CMakeFiles/test_README_example.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/mthind/Desktop/code/micm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../test_README_example"
	cd /Users/mthind/Desktop/code/micm/build/test/tutorial && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_README_example.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/tutorial/CMakeFiles/test_README_example.dir/build: test_README_example
.PHONY : test/tutorial/CMakeFiles/test_README_example.dir/build

test/tutorial/CMakeFiles/test_README_example.dir/clean:
	cd /Users/mthind/Desktop/code/micm/build/test/tutorial && $(CMAKE_COMMAND) -P CMakeFiles/test_README_example.dir/cmake_clean.cmake
.PHONY : test/tutorial/CMakeFiles/test_README_example.dir/clean

test/tutorial/CMakeFiles/test_README_example.dir/depend:
	cd /Users/mthind/Desktop/code/micm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mthind/Desktop/code/micm /Users/mthind/Desktop/code/micm/test/tutorial /Users/mthind/Desktop/code/micm/build /Users/mthind/Desktop/code/micm/build/test/tutorial /Users/mthind/Desktop/code/micm/build/test/tutorial/CMakeFiles/test_README_example.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : test/tutorial/CMakeFiles/test_README_example.dir/depend

