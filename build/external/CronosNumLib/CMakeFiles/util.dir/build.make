# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/dproksch/master_thesis/CronosCelerityBlock

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dproksch/master_thesis/CronosCelerityBlock/build

# Include any dependencies generated for this target.
include external/CronosNumLib/CMakeFiles/util.dir/depend.make

# Include the progress variables for this target.
include external/CronosNumLib/CMakeFiles/util.dir/progress.make

# Include the compile flags for this target's objects.
include external/CronosNumLib/CMakeFiles/util.dir/flags.make

external/CronosNumLib/CMakeFiles/util.dir/util/CronosOstream.C.o: external/CronosNumLib/CMakeFiles/util.dir/flags.make
external/CronosNumLib/CMakeFiles/util.dir/util/CronosOstream.C.o: ../external/CronosNumLib/util/CronosOstream.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dproksch/master_thesis/CronosCelerityBlock/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/CronosNumLib/CMakeFiles/util.dir/util/CronosOstream.C.o"
	cd /home/dproksch/master_thesis/CronosCelerityBlock/build/external/CronosNumLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/util.dir/util/CronosOstream.C.o -c /home/dproksch/master_thesis/CronosCelerityBlock/external/CronosNumLib/util/CronosOstream.C

external/CronosNumLib/CMakeFiles/util.dir/util/CronosOstream.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/util.dir/util/CronosOstream.C.i"
	cd /home/dproksch/master_thesis/CronosCelerityBlock/build/external/CronosNumLib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dproksch/master_thesis/CronosCelerityBlock/external/CronosNumLib/util/CronosOstream.C > CMakeFiles/util.dir/util/CronosOstream.C.i

external/CronosNumLib/CMakeFiles/util.dir/util/CronosOstream.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/util.dir/util/CronosOstream.C.s"
	cd /home/dproksch/master_thesis/CronosCelerityBlock/build/external/CronosNumLib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dproksch/master_thesis/CronosCelerityBlock/external/CronosNumLib/util/CronosOstream.C -o CMakeFiles/util.dir/util/CronosOstream.C.s

external/CronosNumLib/CMakeFiles/util.dir/util/ParameterFileReader.C.o: external/CronosNumLib/CMakeFiles/util.dir/flags.make
external/CronosNumLib/CMakeFiles/util.dir/util/ParameterFileReader.C.o: ../external/CronosNumLib/util/ParameterFileReader.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dproksch/master_thesis/CronosCelerityBlock/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object external/CronosNumLib/CMakeFiles/util.dir/util/ParameterFileReader.C.o"
	cd /home/dproksch/master_thesis/CronosCelerityBlock/build/external/CronosNumLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/util.dir/util/ParameterFileReader.C.o -c /home/dproksch/master_thesis/CronosCelerityBlock/external/CronosNumLib/util/ParameterFileReader.C

external/CronosNumLib/CMakeFiles/util.dir/util/ParameterFileReader.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/util.dir/util/ParameterFileReader.C.i"
	cd /home/dproksch/master_thesis/CronosCelerityBlock/build/external/CronosNumLib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dproksch/master_thesis/CronosCelerityBlock/external/CronosNumLib/util/ParameterFileReader.C > CMakeFiles/util.dir/util/ParameterFileReader.C.i

external/CronosNumLib/CMakeFiles/util.dir/util/ParameterFileReader.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/util.dir/util/ParameterFileReader.C.s"
	cd /home/dproksch/master_thesis/CronosCelerityBlock/build/external/CronosNumLib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dproksch/master_thesis/CronosCelerityBlock/external/CronosNumLib/util/ParameterFileReader.C -o CMakeFiles/util.dir/util/ParameterFileReader.C.s

external/CronosNumLib/CMakeFiles/util.dir/util/util.C.o: external/CronosNumLib/CMakeFiles/util.dir/flags.make
external/CronosNumLib/CMakeFiles/util.dir/util/util.C.o: ../external/CronosNumLib/util/util.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dproksch/master_thesis/CronosCelerityBlock/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object external/CronosNumLib/CMakeFiles/util.dir/util/util.C.o"
	cd /home/dproksch/master_thesis/CronosCelerityBlock/build/external/CronosNumLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/util.dir/util/util.C.o -c /home/dproksch/master_thesis/CronosCelerityBlock/external/CronosNumLib/util/util.C

external/CronosNumLib/CMakeFiles/util.dir/util/util.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/util.dir/util/util.C.i"
	cd /home/dproksch/master_thesis/CronosCelerityBlock/build/external/CronosNumLib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dproksch/master_thesis/CronosCelerityBlock/external/CronosNumLib/util/util.C > CMakeFiles/util.dir/util/util.C.i

external/CronosNumLib/CMakeFiles/util.dir/util/util.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/util.dir/util/util.C.s"
	cd /home/dproksch/master_thesis/CronosCelerityBlock/build/external/CronosNumLib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dproksch/master_thesis/CronosCelerityBlock/external/CronosNumLib/util/util.C -o CMakeFiles/util.dir/util/util.C.s

# Object files for target util
util_OBJECTS = \
"CMakeFiles/util.dir/util/CronosOstream.C.o" \
"CMakeFiles/util.dir/util/ParameterFileReader.C.o" \
"CMakeFiles/util.dir/util/util.C.o"

# External object files for target util
util_EXTERNAL_OBJECTS =

external/CronosNumLib/libutil.a: external/CronosNumLib/CMakeFiles/util.dir/util/CronosOstream.C.o
external/CronosNumLib/libutil.a: external/CronosNumLib/CMakeFiles/util.dir/util/ParameterFileReader.C.o
external/CronosNumLib/libutil.a: external/CronosNumLib/CMakeFiles/util.dir/util/util.C.o
external/CronosNumLib/libutil.a: external/CronosNumLib/CMakeFiles/util.dir/build.make
external/CronosNumLib/libutil.a: external/CronosNumLib/CMakeFiles/util.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/dproksch/master_thesis/CronosCelerityBlock/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library libutil.a"
	cd /home/dproksch/master_thesis/CronosCelerityBlock/build/external/CronosNumLib && $(CMAKE_COMMAND) -P CMakeFiles/util.dir/cmake_clean_target.cmake
	cd /home/dproksch/master_thesis/CronosCelerityBlock/build/external/CronosNumLib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/util.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/CronosNumLib/CMakeFiles/util.dir/build: external/CronosNumLib/libutil.a

.PHONY : external/CronosNumLib/CMakeFiles/util.dir/build

external/CronosNumLib/CMakeFiles/util.dir/clean:
	cd /home/dproksch/master_thesis/CronosCelerityBlock/build/external/CronosNumLib && $(CMAKE_COMMAND) -P CMakeFiles/util.dir/cmake_clean.cmake
.PHONY : external/CronosNumLib/CMakeFiles/util.dir/clean

external/CronosNumLib/CMakeFiles/util.dir/depend:
	cd /home/dproksch/master_thesis/CronosCelerityBlock/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dproksch/master_thesis/CronosCelerityBlock /home/dproksch/master_thesis/CronosCelerityBlock/external/CronosNumLib /home/dproksch/master_thesis/CronosCelerityBlock/build /home/dproksch/master_thesis/CronosCelerityBlock/build/external/CronosNumLib /home/dproksch/master_thesis/CronosCelerityBlock/build/external/CronosNumLib/CMakeFiles/util.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/CronosNumLib/CMakeFiles/util.dir/depend

