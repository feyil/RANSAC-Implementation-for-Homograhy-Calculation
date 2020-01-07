# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/feyil/Desktop/CENG391 Fall 2019/ceng391_06T"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/feyil/Desktop/CENG391 Fall 2019/ceng391_06T/build"

# Utility rule file for test-valgrind.

# Include the progress variables for this target.
include CMakeFiles/test-valgrind.dir/progress.make

CMakeFiles/test-valgrind: image-test
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir="/home/feyil/Desktop/CENG391 Fall 2019/ceng391_06T/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Check for memory leaks using Valgrind"
	valgrind --leak-check=full --show-leak-kinds=all ./image-test

test-valgrind: CMakeFiles/test-valgrind
test-valgrind: CMakeFiles/test-valgrind.dir/build.make

.PHONY : test-valgrind

# Rule to build all files generated by this target.
CMakeFiles/test-valgrind.dir/build: test-valgrind

.PHONY : CMakeFiles/test-valgrind.dir/build

CMakeFiles/test-valgrind.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test-valgrind.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test-valgrind.dir/clean

CMakeFiles/test-valgrind.dir/depend:
	cd "/home/feyil/Desktop/CENG391 Fall 2019/ceng391_06T/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/feyil/Desktop/CENG391 Fall 2019/ceng391_06T" "/home/feyil/Desktop/CENG391 Fall 2019/ceng391_06T" "/home/feyil/Desktop/CENG391 Fall 2019/ceng391_06T/build" "/home/feyil/Desktop/CENG391 Fall 2019/ceng391_06T/build" "/home/feyil/Desktop/CENG391 Fall 2019/ceng391_06T/build/CMakeFiles/test-valgrind.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/test-valgrind.dir/depend

