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
include test/integration/CMakeFiles/test_terminator.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/integration/CMakeFiles/test_terminator.dir/compiler_depend.make

# Include the progress variables for this target.
include test/integration/CMakeFiles/test_terminator.dir/progress.make

# Include the compile flags for this target's objects.
include test/integration/CMakeFiles/test_terminator.dir/flags.make

test/integration/CMakeFiles/test_terminator.dir/test_terminator.cpp.o: test/integration/CMakeFiles/test_terminator.dir/flags.make
test/integration/CMakeFiles/test_terminator.dir/test_terminator.cpp.o: /Users/mthind/Desktop/code/micm/test/integration/test_terminator.cpp
test/integration/CMakeFiles/test_terminator.dir/test_terminator.cpp.o: test/integration/CMakeFiles/test_terminator.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/mthind/Desktop/code/micm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/integration/CMakeFiles/test_terminator.dir/test_terminator.cpp.o"
	cd /Users/mthind/Desktop/code/micm/build/test/integration && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/integration/CMakeFiles/test_terminator.dir/test_terminator.cpp.o -MF CMakeFiles/test_terminator.dir/test_terminator.cpp.o.d -o CMakeFiles/test_terminator.dir/test_terminator.cpp.o -c /Users/mthind/Desktop/code/micm/test/integration/test_terminator.cpp

test/integration/CMakeFiles/test_terminator.dir/test_terminator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test_terminator.dir/test_terminator.cpp.i"
	cd /Users/mthind/Desktop/code/micm/build/test/integration && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mthind/Desktop/code/micm/test/integration/test_terminator.cpp > CMakeFiles/test_terminator.dir/test_terminator.cpp.i

test/integration/CMakeFiles/test_terminator.dir/test_terminator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test_terminator.dir/test_terminator.cpp.s"
	cd /Users/mthind/Desktop/code/micm/build/test/integration && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mthind/Desktop/code/micm/test/integration/test_terminator.cpp -o CMakeFiles/test_terminator.dir/test_terminator.cpp.s

# Object files for target test_terminator
test_terminator_OBJECTS = \
"CMakeFiles/test_terminator.dir/test_terminator.cpp.o"

# External object files for target test_terminator
test_terminator_EXTERNAL_OBJECTS =

test_terminator: test/integration/CMakeFiles/test_terminator.dir/test_terminator.cpp.o
test_terminator: test/integration/CMakeFiles/test_terminator.dir/build.make
test_terminator: lib/libgtest_main.a
test_terminator: _deps/yaml-cpp-build/libyaml-cpp.a
test_terminator: lib/libgtest.a
test_terminator: test/integration/CMakeFiles/test_terminator.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/mthind/Desktop/code/micm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../test_terminator"
	cd /Users/mthind/Desktop/code/micm/build/test/integration && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_terminator.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/integration/CMakeFiles/test_terminator.dir/build: test_terminator
.PHONY : test/integration/CMakeFiles/test_terminator.dir/build

test/integration/CMakeFiles/test_terminator.dir/clean:
	cd /Users/mthind/Desktop/code/micm/build/test/integration && $(CMAKE_COMMAND) -P CMakeFiles/test_terminator.dir/cmake_clean.cmake
.PHONY : test/integration/CMakeFiles/test_terminator.dir/clean

test/integration/CMakeFiles/test_terminator.dir/depend:
	cd /Users/mthind/Desktop/code/micm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mthind/Desktop/code/micm /Users/mthind/Desktop/code/micm/test/integration /Users/mthind/Desktop/code/micm/build /Users/mthind/Desktop/code/micm/build/test/integration /Users/mthind/Desktop/code/micm/build/test/integration/CMakeFiles/test_terminator.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : test/integration/CMakeFiles/test_terminator.dir/depend

