# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/local/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/andi/SW/libpointmesher

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/andi/SW/libpointmesher/build

# Include any dependencies generated for this target.
include CMakeFiles/meshFilter.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/meshFilter.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/meshFilter.dir/flags.make

CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o: CMakeFiles/meshFilter.dir/flags.make
CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o: ../tests/meshFilter.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/andi/SW/libpointmesher/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o"
	/home/andi/sbin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o -c /home/andi/SW/libpointmesher/tests/meshFilter.cpp

CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.i"
	/home/andi/sbin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/andi/SW/libpointmesher/tests/meshFilter.cpp > CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.i

CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.s"
	/home/andi/sbin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/andi/SW/libpointmesher/tests/meshFilter.cpp -o CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.s

CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o.requires:
.PHONY : CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o.requires

CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o.provides: CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o.requires
	$(MAKE) -f CMakeFiles/meshFilter.dir/build.make CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o.provides.build
.PHONY : CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o.provides

CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o.provides.build: CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o

# Object files for target meshFilter
meshFilter_OBJECTS = \
"CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o"

# External object files for target meshFilter
meshFilter_EXTERNAL_OBJECTS =

meshFilter: CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o
meshFilter: libpointmesher.a
meshFilter: /home/andi/SW/libpointmatcher/build/libpointmatcher.a
meshFilter: /home/andi/SW/libs/libnabo/build/libnabo.a
meshFilter: /home/andi/SW/OpenMesh-2.0/build/Build/lib/OpenMesh/libOpenMeshCore.so
meshFilter: /home/andi/SW/IsoEx/build/libisoex.a
meshFilter: /home/andi/SW/libs/flann-1.6.8-src/lib/libflann_cpp_s.a
meshFilter: /home/andi/SW/libs/qhull-2010.1/src2/qhull/libqhull.a
meshFilter: /usr/local/lib/libpcl_io.so
meshFilter: /usr/local/lib/libpcl_filters.so
meshFilter: /usr/local/lib/libpcl_features.so
meshFilter: /usr/local/lib/libpcl_segmentation.so
meshFilter: /usr/local/lib/libpcl_surface.so
meshFilter: /usr/local/lib/libpcl_registration.so
meshFilter: /usr/local/lib/libCGAL.so
meshFilter: /usr/lib/libOpenCL.so
meshFilter: CMakeFiles/meshFilter.dir/build.make
meshFilter: CMakeFiles/meshFilter.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable meshFilter"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/meshFilter.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/meshFilter.dir/build: meshFilter
.PHONY : CMakeFiles/meshFilter.dir/build

CMakeFiles/meshFilter.dir/requires: CMakeFiles/meshFilter.dir/tests/meshFilter.cpp.o.requires
.PHONY : CMakeFiles/meshFilter.dir/requires

CMakeFiles/meshFilter.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/meshFilter.dir/cmake_clean.cmake
.PHONY : CMakeFiles/meshFilter.dir/clean

CMakeFiles/meshFilter.dir/depend:
	cd /home/andi/SW/libpointmesher/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andi/SW/libpointmesher /home/andi/SW/libpointmesher /home/andi/SW/libpointmesher/build /home/andi/SW/libpointmesher/build /home/andi/SW/libpointmesher/build/CMakeFiles/meshFilter.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/meshFilter.dir/depend
