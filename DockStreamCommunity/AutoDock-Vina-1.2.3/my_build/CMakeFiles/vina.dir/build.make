# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

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
$$(VERBOSE).SILENT:


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
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/home/linjie/CLionProjects/AutoDock-Vina-1.2.3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/home/linjie/CLionProjects/AutoDock-Vina-1.2.3/my_build

# Include any dependencies generated for this target.
include CMakeFiles/vina.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/vina.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/vina.dir/flags.make

CMakeFiles/vina.dir/src/main/main.cpp.o: CMakeFiles/vina.dir/flags.make
CMakeFiles/vina.dir/src/main/main.cpp.o: ../src/main/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/home/linjie/CLionProjects/AutoDock-Vina-1.2.3/my_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/vina.dir/src/main/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vina.dir/src/main/main.cpp.o -c /mnt/home/linjie/CLionProjects/AutoDock-Vina-1.2.3/src/main/main.cpp

CMakeFiles/vina.dir/src/main/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vina.dir/src/main/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/home/linjie/CLionProjects/AutoDock-Vina-1.2.3/src/main/main.cpp > CMakeFiles/vina.dir/src/main/main.cpp.i

CMakeFiles/vina.dir/src/main/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vina.dir/src/main/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/home/linjie/CLionProjects/AutoDock-Vina-1.2.3/src/main/main.cpp -o CMakeFiles/vina.dir/src/main/main.cpp.s

# Object files for target vina
vina_OBJECTS = \
"CMakeFiles/vina.dir/src/main/main.cpp.o"

# External object files for target vina
vina_EXTERNAL_OBJECTS =

vina: CMakeFiles/vina.dir/src/main/main.cpp.o
vina: CMakeFiles/vina.dir/build.make
vina: libvina_lib.a
vina: CMakeFiles/vina.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/home/linjie/CLionProjects/AutoDock-Vina-1.2.3/my_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable vina"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vina.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/vina.dir/build: vina

.PHONY : CMakeFiles/vina.dir/build

CMakeFiles/vina.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/vina.dir/cmake_clean.cmake
.PHONY : CMakeFiles/vina.dir/clean

CMakeFiles/vina.dir/depend:
	cd /mnt/home/linjie/CLionProjects/AutoDock-Vina-1.2.3/my_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/home/linjie/CLionProjects/AutoDock-Vina-1.2.3 /mnt/home/linjie/CLionProjects/AutoDock-Vina-1.2.3 /mnt/home/linjie/CLionProjects/AutoDock-Vina-1.2.3/my_build /mnt/home/linjie/CLionProjects/AutoDock-Vina-1.2.3/my_build /mnt/home/linjie/CLionProjects/AutoDock-Vina-1.2.3/my_build/CMakeFiles/vina.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/vina.dir/depend
