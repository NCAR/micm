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
include test/unit/solver/CMakeFiles/test_solver_builder.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/unit/solver/CMakeFiles/test_solver_builder.dir/compiler_depend.make

# Include the progress variables for this target.
include test/unit/solver/CMakeFiles/test_solver_builder.dir/progress.make

# Include the compile flags for this target's objects.
include test/unit/solver/CMakeFiles/test_solver_builder.dir/flags.make

test/unit/solver/CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.o: test/unit/solver/CMakeFiles/test_solver_builder.dir/flags.make
test/unit/solver/CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.o: /Users/mthind/Desktop/code/micm/test/unit/solver/test_solver_builder.cpp
test/unit/solver/CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.o: test/unit/solver/CMakeFiles/test_solver_builder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/mthind/Desktop/code/micm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/unit/solver/CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.o"
	cd /Users/mthind/Desktop/code/micm/build/test/unit/solver && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/unit/solver/CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.o -MF CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.o.d -o CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.o -c /Users/mthind/Desktop/code/micm/test/unit/solver/test_solver_builder.cpp

test/unit/solver/CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.i"
	cd /Users/mthind/Desktop/code/micm/build/test/unit/solver && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mthind/Desktop/code/micm/test/unit/solver/test_solver_builder.cpp > CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.i

test/unit/solver/CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.s"
	cd /Users/mthind/Desktop/code/micm/build/test/unit/solver && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mthind/Desktop/code/micm/test/unit/solver/test_solver_builder.cpp -o CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.s

# Object files for target test_solver_builder
test_solver_builder_OBJECTS = \
"CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.o"

# External object files for target test_solver_builder
test_solver_builder_EXTERNAL_OBJECTS =

test_solver_builder: test/unit/solver/CMakeFiles/test_solver_builder.dir/test_solver_builder.cpp.o
test_solver_builder: test/unit/solver/CMakeFiles/test_solver_builder.dir/build.make
test_solver_builder: lib/libgtest_main.a
test_solver_builder: _deps/yaml-cpp-build/libyaml-cpp.a
test_solver_builder: lib/libgtest.a
test_solver_builder: test/unit/solver/CMakeFiles/test_solver_builder.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/mthind/Desktop/code/micm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../test_solver_builder"
	cd /Users/mthind/Desktop/code/micm/build/test/unit/solver && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_solver_builder.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/unit/solver/CMakeFiles/test_solver_builder.dir/build: test_solver_builder
.PHONY : test/unit/solver/CMakeFiles/test_solver_builder.dir/build

test/unit/solver/CMakeFiles/test_solver_builder.dir/clean:
	cd /Users/mthind/Desktop/code/micm/build/test/unit/solver && $(CMAKE_COMMAND) -P CMakeFiles/test_solver_builder.dir/cmake_clean.cmake
.PHONY : test/unit/solver/CMakeFiles/test_solver_builder.dir/clean

test/unit/solver/CMakeFiles/test_solver_builder.dir/depend:
	cd /Users/mthind/Desktop/code/micm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mthind/Desktop/code/micm /Users/mthind/Desktop/code/micm/test/unit/solver /Users/mthind/Desktop/code/micm/build /Users/mthind/Desktop/code/micm/build/test/unit/solver /Users/mthind/Desktop/code/micm/build/test/unit/solver/CMakeFiles/test_solver_builder.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : test/unit/solver/CMakeFiles/test_solver_builder.dir/depend

