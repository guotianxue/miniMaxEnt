# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

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
CMAKE_COMMAND = D:/cmake/bin/cmake.exe

# The command to remove a file.
RM = D:/cmake/bin/cmake.exe -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = E:/leilei/code_use

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = E:/leilei/code_use/build

# Include any dependencies generated for this target.
include CMakeFiles/Raw_data_Inverse.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/Raw_data_Inverse.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Raw_data_Inverse.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Raw_data_Inverse.dir/flags.make

CMakeFiles/Raw_data_Inverse.dir/test.cpp.obj: CMakeFiles/Raw_data_Inverse.dir/flags.make
CMakeFiles/Raw_data_Inverse.dir/test.cpp.obj: CMakeFiles/Raw_data_Inverse.dir/includes_CXX.rsp
CMakeFiles/Raw_data_Inverse.dir/test.cpp.obj: ../test.cpp
CMakeFiles/Raw_data_Inverse.dir/test.cpp.obj: CMakeFiles/Raw_data_Inverse.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=E:/leilei/code_use/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Raw_data_Inverse.dir/test.cpp.obj"
	D:/MinGW/bin/g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Raw_data_Inverse.dir/test.cpp.obj -MF CMakeFiles/Raw_data_Inverse.dir/test.cpp.obj.d -o CMakeFiles/Raw_data_Inverse.dir/test.cpp.obj -c E:/leilei/code_use/test.cpp

CMakeFiles/Raw_data_Inverse.dir/test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Raw_data_Inverse.dir/test.cpp.i"
	D:/MinGW/bin/g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E E:/leilei/code_use/test.cpp > CMakeFiles/Raw_data_Inverse.dir/test.cpp.i

CMakeFiles/Raw_data_Inverse.dir/test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Raw_data_Inverse.dir/test.cpp.s"
	D:/MinGW/bin/g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S E:/leilei/code_use/test.cpp -o CMakeFiles/Raw_data_Inverse.dir/test.cpp.s

# Object files for target Raw_data_Inverse
Raw_data_Inverse_OBJECTS = \
"CMakeFiles/Raw_data_Inverse.dir/test.cpp.obj"

# External object files for target Raw_data_Inverse
Raw_data_Inverse_EXTERNAL_OBJECTS =

Raw_data_Inverse.exe: CMakeFiles/Raw_data_Inverse.dir/test.cpp.obj
Raw_data_Inverse.exe: CMakeFiles/Raw_data_Inverse.dir/build.make
Raw_data_Inverse.exe: CMakeFiles/Raw_data_Inverse.dir/linklibs.rsp
Raw_data_Inverse.exe: CMakeFiles/Raw_data_Inverse.dir/objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=E:/leilei/code_use/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Raw_data_Inverse.exe"
	D:/cmake/bin/cmake.exe -E rm -f CMakeFiles/Raw_data_Inverse.dir/objects.a
	D:/MinGW/bin/ar.exe qc CMakeFiles/Raw_data_Inverse.dir/objects.a @CMakeFiles/Raw_data_Inverse.dir/objects1.rsp
	D:/MinGW/bin/g++.exe -g -Wl,--whole-archive CMakeFiles/Raw_data_Inverse.dir/objects.a -Wl,--no-whole-archive -o Raw_data_Inverse.exe -Wl,--out-implib,libRaw_data_Inverse.dll.a -Wl,--major-image-version,0,--minor-image-version,0 @CMakeFiles/Raw_data_Inverse.dir/linklibs.rsp

# Rule to build all files generated by this target.
CMakeFiles/Raw_data_Inverse.dir/build: Raw_data_Inverse.exe
.PHONY : CMakeFiles/Raw_data_Inverse.dir/build

CMakeFiles/Raw_data_Inverse.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Raw_data_Inverse.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Raw_data_Inverse.dir/clean

CMakeFiles/Raw_data_Inverse.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" E:/leilei/code_use E:/leilei/code_use E:/leilei/code_use/build E:/leilei/code_use/build E:/leilei/code_use/build/CMakeFiles/Raw_data_Inverse.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Raw_data_Inverse.dir/depend
