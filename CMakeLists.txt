# Declare minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

project(sod2d)           # Project name
enable_language(Fortran) # Supported language Fortran

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
	"Choose the installation directory; by default it installs in the SOD2D directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-fopenmp -cpp -funroll-all-loops -fno-f2c -O3 -ffree-line-length-512 -march=native")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fopenmp -cpp -ffree-line-length-512 -fbounds-check -fbacktrace -ftrapv -Wconversion-extra -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-mtune=native -list-line-len=512 -fpp -mtune=native -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-mtune=native -list-line-len=512 -fpp -ftrapuv -check all,noarg_temp_created -traceback -debug full -warn all,nodec,nointerfaces -fp-stack-check -ansi-alias -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "nvfortran")
  # nvfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-fast -D_USE_NVTX -cpp -gpu=cc70,managed,lineinfo -cuda -acc -Minfo=accel")
  #set (CMAKE_Fortran_FLAGS_RELEASE "-fast -D_USE_NVTX -cpp -cuda -Minfo=accel")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -DGPU -cpp -Mcuda=cc70 -Wall -C -Ktrap=fp -g")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

enable_testing()
ADD_SUBDIRECTORY(unitt)

# Compile
file(GLOB sources src/*.f90 src/main/*.f90) # Adds all files to compilation list, with dependencies pre analyzed by CMake
add_executable(sod2d.x ${sources})   # Build executable from listed sources
