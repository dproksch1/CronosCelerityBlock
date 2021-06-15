# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity

# Include any dependencies generated for this target.
include external/CronosNumLib/CMakeFiles/matrix.dir/depend.make

# Include the progress variables for this target.
include external/CronosNumLib/CMakeFiles/matrix.dir/progress.make

# Include the compile flags for this target's objects.
include external/CronosNumLib/CMakeFiles/matrix.dir/flags.make

external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/array.C.o: external/CronosNumLib/CMakeFiles/matrix.dir/flags.make
external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/array.C.o: external/CronosNumLib/Matrix/array.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/array.C.o"
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/Matrix/array.C.o -c /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib/Matrix/array.C

external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/array.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/Matrix/array.C.i"
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib/Matrix/array.C > CMakeFiles/matrix.dir/Matrix/array.C.i

external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/array.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/Matrix/array.C.s"
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib/Matrix/array.C -o CMakeFiles/matrix.dir/Matrix/array.C.s

external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/boundary.C.o: external/CronosNumLib/CMakeFiles/matrix.dir/flags.make
external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/boundary.C.o: external/CronosNumLib/Matrix/boundary.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/boundary.C.o"
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/Matrix/boundary.C.o -c /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib/Matrix/boundary.C

external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/boundary.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/Matrix/boundary.C.i"
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib/Matrix/boundary.C > CMakeFiles/matrix.dir/Matrix/boundary.C.i

external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/boundary.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/Matrix/boundary.C.s"
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib/Matrix/boundary.C -o CMakeFiles/matrix.dir/Matrix/boundary.C.s

external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/fftmatrix.C.o: external/CronosNumLib/CMakeFiles/matrix.dir/flags.make
external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/fftmatrix.C.o: external/CronosNumLib/Matrix/fftmatrix.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/fftmatrix.C.o"
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/Matrix/fftmatrix.C.o -c /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib/Matrix/fftmatrix.C

external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/fftmatrix.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/Matrix/fftmatrix.C.i"
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib/Matrix/fftmatrix.C > CMakeFiles/matrix.dir/Matrix/fftmatrix.C.i

external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/fftmatrix.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/Matrix/fftmatrix.C.s"
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib/Matrix/fftmatrix.C -o CMakeFiles/matrix.dir/Matrix/fftmatrix.C.s

external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/matrix.C.o: external/CronosNumLib/CMakeFiles/matrix.dir/flags.make
external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/matrix.C.o: external/CronosNumLib/Matrix/matrix.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/matrix.C.o"
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/Matrix/matrix.C.o -c /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib/Matrix/matrix.C

external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/matrix.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/Matrix/matrix.C.i"
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib/Matrix/matrix.C > CMakeFiles/matrix.dir/Matrix/matrix.C.i

external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/matrix.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/Matrix/matrix.C.s"
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib/Matrix/matrix.C -o CMakeFiles/matrix.dir/Matrix/matrix.C.s

# Object files for target matrix
matrix_OBJECTS = \
"CMakeFiles/matrix.dir/Matrix/array.C.o" \
"CMakeFiles/matrix.dir/Matrix/boundary.C.o" \
"CMakeFiles/matrix.dir/Matrix/fftmatrix.C.o" \
"CMakeFiles/matrix.dir/Matrix/matrix.C.o"

# External object files for target matrix
matrix_EXTERNAL_OBJECTS =

external/CronosNumLib/libmatrix.a: external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/array.C.o
external/CronosNumLib/libmatrix.a: external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/boundary.C.o
external/CronosNumLib/libmatrix.a: external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/fftmatrix.C.o
external/CronosNumLib/libmatrix.a: external/CronosNumLib/CMakeFiles/matrix.dir/Matrix/matrix.C.o
external/CronosNumLib/libmatrix.a: external/CronosNumLib/CMakeFiles/matrix.dir/build.make
external/CronosNumLib/libmatrix.a: external/CronosNumLib/CMakeFiles/matrix.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX static library libmatrix.a"
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && $(CMAKE_COMMAND) -P CMakeFiles/matrix.dir/cmake_clean_target.cmake
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/matrix.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/CronosNumLib/CMakeFiles/matrix.dir/build: external/CronosNumLib/libmatrix.a

.PHONY : external/CronosNumLib/CMakeFiles/matrix.dir/build

external/CronosNumLib/CMakeFiles/matrix.dir/clean:
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib && $(CMAKE_COMMAND) -P CMakeFiles/matrix.dir/cmake_clean.cmake
.PHONY : external/CronosNumLib/CMakeFiles/matrix.dir/clean

external/CronosNumLib/CMakeFiles/matrix.dir/depend:
	cd /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity/external/CronosNumLib/CMakeFiles/matrix.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/CronosNumLib/CMakeFiles/matrix.dir/depend

