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
include test/unit/system/CMakeFiles/test_phase.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/unit/system/CMakeFiles/test_phase.dir/compiler_depend.make

# Include the progress variables for this target.
include test/unit/system/CMakeFiles/test_phase.dir/progress.make

# Include the compile flags for this target's objects.
include test/unit/system/CMakeFiles/test_phase.dir/flags.make

test/unit/system/CMakeFiles/test_phase.dir/test_phase.cpp.o: test/unit/system/CMakeFiles/test_phase.dir/flags.make
test/unit/system/CMakeFiles/test_phase.dir/test_phase.cpp.o: /Users/mthind/Desktop/code/micm/test/unit/system/test_phase.cpp
test/unit/system/CMakeFiles/test_phase.dir/test_phase.cpp.o: test/unit/system/CMakeFiles/test_phase.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/mthind/Desktop/code/micm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/unit/system/CMakeFiles/test_phase.dir/test_phase.cpp.o"
	cd /Users/mthind/Desktop/code/micm/build/test/unit/system && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/unit/system/CMakeFiles/test_phase.dir/test_phase.cpp.o -MF CMakeFiles/test_phase.dir/test_phase.cpp.o.d -o CMakeFiles/test_phase.dir/test_phase.cpp.o -c /Users/mthind/Desktop/code/micm/test/unit/system/test_phase.cpp

test/unit/system/CMakeFiles/test_phase.dir/test_phase.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test_phase.dir/test_phase.cpp.i"
	cd /Users/mthind/Desktop/code/micm/build/test/unit/system && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mthind/Desktop/code/micm/test/unit/system/test_phase.cpp > CMakeFiles/test_phase.dir/test_phase.cpp.i

test/unit/system/CMakeFiles/test_phase.dir/test_phase.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test_phase.dir/test_phase.cpp.s"
	cd /Users/mthind/Desktop/code/micm/build/test/unit/system && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mthind/Desktop/code/micm/test/unit/system/test_phase.cpp -o CMakeFiles/test_phase.dir/test_phase.cpp.s

# Object files for target test_phase
test_phase_OBJECTS = \
"CMakeFiles/test_phase.dir/test_phase.cpp.o"

# External object files for target test_phase
test_phase_EXTERNAL_OBJECTS =

test_phase: test/unit/system/CMakeFiles/test_phase.dir/test_phase.cpp.o
test_phase: test/unit/system/CMakeFiles/test_phase.dir/build.make
test_phase: lib/libgtest_main.a
test_phase: _deps/yaml-cpp-build/libyaml-cpp.a
test_phase: lib/libgtest.a
test_phase: test/unit/system/CMakeFiles/test_phase.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/mthind/Desktop/code/micm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../test_phase"
	cd /Users/mthind/Desktop/code/micm/build/test/unit/system && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_phase.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/unit/system/CMakeFiles/test_phase.dir/build: test_phase
.PHONY : test/unit/system/CMakeFiles/test_phase.dir/build

test/unit/system/CMakeFiles/test_phase.dir/clean:
	cd /Users/mthind/Desktop/code/micm/build/test/unit/system && $(CMAKE_COMMAND) -P CMakeFiles/test_phase.dir/cmake_clean.cmake
.PHONY : test/unit/system/CMakeFiles/test_phase.dir/clean

test/unit/system/CMakeFiles/test_phase.dir/depend:
	cd /Users/mthind/Desktop/code/micm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mthind/Desktop/code/micm /Users/mthind/Desktop/code/micm/test/unit/system /Users/mthind/Desktop/code/micm/build /Users/mthind/Desktop/code/micm/build/test/unit/system /Users/mthind/Desktop/code/micm/build/test/unit/system/CMakeFiles/test_phase.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : test/unit/system/CMakeFiles/test_phase.dir/depend

