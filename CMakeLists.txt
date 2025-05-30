# This is the CMake build file for the MPASSIT
#
# Larissa Reames

cmake_minimum_required(VERSION 3.15)

# Get the version from the VERSION file.
file(STRINGS "VERSION" pVersion)

project(mpassit
  VERSION ${pVersion}
  LANGUAGES C Fortran )
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# User options.
option(OPENMP "use OpenMP threading" ON)
option(ENABLE_DOCS "Enable generation of doxygen-based documentation." OFF)

# Set compiler flags.
if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)(LLVM)?$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model precise")
  if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fp-stack-check -fstack-protector-all -fpe0 -debug -ftrapuv")
  else()
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check bounds -check contiguous -check format -check stack -check shape -check assume -fstack-protector-all -fpe0 -debug -ftrapuv")
  endif()
  if(APPLE)
    # The linker on macOS does not include `common symbols` by default.
    # Passing the -c flag includes them and fixes an error with undefined symbols.
    set(CMAKE_Fortran_ARCHIVE_FINISH "<CMAKE_RANLIB> -c <TARGET>")
    set(CMAKE_C_ARCHIVE_FINISH "<CMAKE_RANLIB> -c <TARGET>")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fbacktrace")
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fallow-invalid-boz")
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls")
endif()

if(CMAKE_C_COMPILER_ID MATCHES "^(Intel)(LLVM)?$")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -traceback")
  set(CMAKE_C_FLAGS_RELEASE "-O2")
  set(CMAKE_C_FLAGS_DEBUG "-O0")
endif()

# Find packages.
find_package(NetCDF 4.3.3 REQUIRED C Fortran)
find_package(MPI REQUIRED C Fortran)
find_package(ESMF 8.3.0 REQUIRED)

if(OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

# If doxygen documentation we enabled, build it. This must come before
# adding the source code directories; the main documentation build
# must happen before any of the utility document builds.
if(ENABLE_DOCS)
  find_package(Doxygen REQUIRED)
  set(abs_top_srcdir "${CMAKE_SOURCE_DIR}")
  add_subdirectory(docs)  
endif()

set(lib_src
    model_grid.F90
    program_setup.F90
    utils.F90
    input_data.F90
    interp.F90
    write_data.F90
    datetime_module.F90
    constants_module.F90
    llxy_module.F90
    misc_definitions_module.F90
    module_map_utils.F90
)

set(exe_src mpassit.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)(LLVM)?$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
  
  # Turn on this argument mismatch flag for gfortran10.
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()
endif()

set(exe_name mpassit)

add_library(mpassit_lib STATIC ${lib_src})
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(mpassit_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(mpassit_lib INTERFACE ${mod_dir})

target_link_libraries(
  mpassit_lib
  PUBLIC
  esmf
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(${exe_name} PUBLIC OpenMP::OpenMP_Fortran)
endif()

target_link_libraries(${exe_name} PRIVATE mpassit_lib)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()

