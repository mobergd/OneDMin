# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /blues/gpfs/home/kmoore/miniconda/envs/pacc-env-mini/bin/cmake

# The command to remove a file.
RM = /blues/gpfs/home/kmoore/miniconda/envs/pacc-env-mini/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /lcrc/project/CMRP/pacc/OneDMin

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /lcrc/project/CMRP/pacc/OneDMin/build

# Include any dependencies generated for this target.
include CMakeFiles/ljlibs.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ljlibs.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ljlibs.dir/flags.make

CMakeFiles/ljlibs.dir/src/mass.f.o: CMakeFiles/ljlibs.dir/flags.make
CMakeFiles/ljlibs.dir/src/mass.f.o: ../src/mass.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lcrc/project/CMRP/pacc/OneDMin/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/ljlibs.dir/src/mass.f.o"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /lcrc/project/CMRP/pacc/OneDMin/src/mass.f -o CMakeFiles/ljlibs.dir/src/mass.f.o

CMakeFiles/ljlibs.dir/src/mass.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/ljlibs.dir/src/mass.f.i"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /lcrc/project/CMRP/pacc/OneDMin/src/mass.f > CMakeFiles/ljlibs.dir/src/mass.f.i

CMakeFiles/ljlibs.dir/src/mass.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/ljlibs.dir/src/mass.f.s"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /lcrc/project/CMRP/pacc/OneDMin/src/mass.f -o CMakeFiles/ljlibs.dir/src/mass.f.s

CMakeFiles/ljlibs.dir/src/spin.f.o: CMakeFiles/ljlibs.dir/flags.make
CMakeFiles/ljlibs.dir/src/spin.f.o: ../src/spin.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lcrc/project/CMRP/pacc/OneDMin/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object CMakeFiles/ljlibs.dir/src/spin.f.o"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /lcrc/project/CMRP/pacc/OneDMin/src/spin.f -o CMakeFiles/ljlibs.dir/src/spin.f.o

CMakeFiles/ljlibs.dir/src/spin.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/ljlibs.dir/src/spin.f.i"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /lcrc/project/CMRP/pacc/OneDMin/src/spin.f > CMakeFiles/ljlibs.dir/src/spin.f.i

CMakeFiles/ljlibs.dir/src/spin.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/ljlibs.dir/src/spin.f.s"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /lcrc/project/CMRP/pacc/OneDMin/src/spin.f -o CMakeFiles/ljlibs.dir/src/spin.f.s

# Object files for target ljlibs
ljlibs_OBJECTS = \
"CMakeFiles/ljlibs.dir/src/mass.f.o" \
"CMakeFiles/ljlibs.dir/src/spin.f.o"

# External object files for target ljlibs
ljlibs_EXTERNAL_OBJECTS =

libljlibs.a: CMakeFiles/ljlibs.dir/src/mass.f.o
libljlibs.a: CMakeFiles/ljlibs.dir/src/spin.f.o
libljlibs.a: CMakeFiles/ljlibs.dir/build.make
libljlibs.a: CMakeFiles/ljlibs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/lcrc/project/CMRP/pacc/OneDMin/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking Fortran static library libljlibs.a"
	$(CMAKE_COMMAND) -P CMakeFiles/ljlibs.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ljlibs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ljlibs.dir/build: libljlibs.a

.PHONY : CMakeFiles/ljlibs.dir/build

CMakeFiles/ljlibs.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ljlibs.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ljlibs.dir/clean

CMakeFiles/ljlibs.dir/depend:
	cd /lcrc/project/CMRP/pacc/OneDMin/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lcrc/project/CMRP/pacc/OneDMin /lcrc/project/CMRP/pacc/OneDMin /lcrc/project/CMRP/pacc/OneDMin/build /lcrc/project/CMRP/pacc/OneDMin/build /lcrc/project/CMRP/pacc/OneDMin/build/CMakeFiles/ljlibs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ljlibs.dir/depend

