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
include test/unit/util/CMakeFiles/test_vector_matrix.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/unit/util/CMakeFiles/test_vector_matrix.dir/compiler_depend.make

# Include the progress variables for this target.
include test/unit/util/CMakeFiles/test_vector_matrix.dir/progress.make

# Include the compile flags for this target's objects.
include test/unit/util/CMakeFiles/test_vector_matrix.dir/flags.make

test/unit/util/CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.o: test/unit/util/CMakeFiles/test_vector_matrix.dir/flags.make
test/unit/util/CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.o: /Users/mthind/Desktop/code/micm/test/unit/util/test_vector_matrix.cpp
test/unit/util/CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.o: test/unit/util/CMakeFiles/test_vector_matrix.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/mthind/Desktop/code/micm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/unit/util/CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.o"
	cd /Users/mthind/Desktop/code/micm/build/test/unit/util && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/unit/util/CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.o -MF CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.o.d -o CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.o -c /Users/mthind/Desktop/code/micm/test/unit/util/test_vector_matrix.cpp

test/unit/util/CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.i"
	cd /Users/mthind/Desktop/code/micm/build/test/unit/util && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mthind/Desktop/code/micm/test/unit/util/test_vector_matrix.cpp > CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.i

test/unit/util/CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.s"
	cd /Users/mthind/Desktop/code/micm/build/test/unit/util && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mthind/Desktop/code/micm/test/unit/util/test_vector_matrix.cpp -o CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.s

# Object files for target test_vector_matrix
test_vector_matrix_OBJECTS = \
"CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.o"

# External object files for target test_vector_matrix
test_vector_matrix_EXTERNAL_OBJECTS =

test_vector_matrix: test/unit/util/CMakeFiles/test_vector_matrix.dir/test_vector_matrix.cpp.o
test_vector_matrix: test/unit/util/CMakeFiles/test_vector_matrix.dir/build.make
test_vector_matrix: lib/libgtest_main.a
test_vector_matrix: _deps/yaml-cpp-build/libyaml-cpp.a
test_vector_matrix: lib/libgtest.a
test_vector_matrix: test/unit/util/CMakeFiles/test_vector_matrix.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/mthind/Desktop/code/micm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../test_vector_matrix"
	cd /Users/mthind/Desktop/code/micm/build/test/unit/util && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_vector_matrix.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/unit/util/CMakeFiles/test_vector_matrix.dir/build: test_vector_matrix
.PHONY : test/unit/util/CMakeFiles/test_vector_matrix.dir/build

test/unit/util/CMakeFiles/test_vector_matrix.dir/clean:
	cd /Users/mthind/Desktop/code/micm/build/test/unit/util && $(CMAKE_COMMAND) -P CMakeFiles/test_vector_matrix.dir/cmake_clean.cmake
.PHONY : test/unit/util/CMakeFiles/test_vector_matrix.dir/clean

test/unit/util/CMakeFiles/test_vector_matrix.dir/depend:
	cd /Users/mthind/Desktop/code/micm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mthind/Desktop/code/micm /Users/mthind/Desktop/code/micm/test/unit/util /Users/mthind/Desktop/code/micm/build /Users/mthind/Desktop/code/micm/build/test/unit/util /Users/mthind/Desktop/code/micm/build/test/unit/util/CMakeFiles/test_vector_matrix.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : test/unit/util/CMakeFiles/test_vector_matrix.dir/depend

