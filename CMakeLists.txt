# Declare minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

project(sod2d)           # Project name
enable_language(Fortran) # Supported language Fortran

# Set GNU Gfortran compiler
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	# Set some compilation flags to variables
	set(dialect "-O3 -march=native -ffree-form -std=f2018 -fimplicit-none")
	set(bounds "-fbounds-check")
endif()

# Set PGI compiler
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    # Set some compilation flags to variables
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}") # Debug exclusive options
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")            # Used for both Release and Debug

#
# Compile
#
file(GLOB_RECURSE sources src/*.f90) # Adds all files to compilation list, with dependencies pre analyzed by CMake
add_executable(sod2d.x ${sources})      # Build executable from listed sources
