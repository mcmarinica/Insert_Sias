# Build of Insert_Sias with CMake.
#
#   cmake -B build
#   cmake --build build
#
# The executables are written to Src/bin, where the examples and the
# documentation expect them. See doc/install.md.

cmake_minimum_required(VERSION 3.16)

# Compiler preference: gfortran first, then ifx, then the classic ifort.
# An explicit choice always wins:
#   cmake -B build -DCMAKE_Fortran_COMPILER=ifx      or      FC=ifx cmake -B build
if(NOT CMAKE_Fortran_COMPILER AND NOT DEFINED ENV{FC})
  find_program(INSERT_SIAS_PREFERRED_FC NAMES gfortran ifx ifort)
  mark_as_advanced(INSERT_SIAS_PREFERRED_FC)
  if(INSERT_SIAS_PREFERRED_FC)
    set(CMAKE_Fortran_COMPILER "${INSERT_SIAS_PREFERRED_FC}"
        CACHE FILEPATH "Fortran compiler")
  endif()
endif()

project(Insert_Sias LANGUAGES Fortran)

# Optimised build unless another type is requested (-DCMAKE_BUILD_TYPE=Debug).
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug)
endif()

# Where the executables go. The default matches the makefiles.
set(INSERT_SIAS_BIN_DIR "${PROJECT_SOURCE_DIR}/Src/bin"
    CACHE PATH "Directory that receives the executables")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${INSERT_SIAS_BIN_DIR}")
set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")

# Compiler-specific flags.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # Some free-form lines exceed 132 columns, the limit of the standard before
  # Fortran 2023. Accept them whatever the default of the gfortran release.
  add_compile_options(-ffree-line-length-none)
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fbacktrace -fcheck=bounds")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")   # ifx (IntelLLVM) and ifort (Intel)
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback -check bounds")
endif()

message(STATUS "Insert_Sias: Fortran compiler ${CMAKE_Fortran_COMPILER_ID} "
               "${CMAKE_Fortran_COMPILER_VERSION} (${CMAKE_Fortran_COMPILER})")
message(STATUS "Insert_Sias: build type ${CMAKE_BUILD_TYPE}, "
               "executables in ${INSERT_SIAS_BIN_DIR}")

set(SRC "${PROJECT_SOURCE_DIR}/Src")

# The routines of the defect builder, Src/lib. Modules come first for the reader;
# CMake works out the compilation order by itself.
add_library(insert_sias STATIC
  ${SRC}/lib/module_defect_types.f90
  ${SRC}/lib/math.f90
  ${SRC}/lib/module_pbc.f90
  ${SRC}/lib/module_octahedre.f90
  ${SRC}/lib/module_migration_sia.f90
  ${SRC}/lib/module_111_smooth.f90
  ${SRC}/lib/module_100_smooth.f90
  ${SRC}/lib/lammps_output.f90
  ${SRC}/lib/input.f90
  ${SRC}/lib/sites.f90
  ${SRC}/lib/point_defects.f90
  ${SRC}/lib/migration.f90
  ${SRC}/lib/clusters.f90
  ${SRC}/lib/kazuto_style.f90
  ${SRC}/lib/rotations.f90
  ${SRC}/lib/indexx.f90
  ${SRC}/lib/lapack.f
)

# One executable per program of Src/programs, named <name>.x as with the makefiles.
function(insert_sias_program name)
  add_executable(${name} ${ARGN})
  set_target_properties(${name} PROPERTIES SUFFIX ".x")
  install(TARGETS ${name} RUNTIME DESTINATION bin)
endfunction()

insert_sias_program(put_vct        ${SRC}/programs/put_vct.f90)
target_link_libraries(put_vct PRIVATE insert_sias)

insert_sias_program(pre_put_vct    ${SRC}/programs/pre_put_vct.f90)
insert_sias_program(cube           ${SRC}/programs/cube.f)
insert_sias_program(cube_shell     ${SRC}/programs/cube_shell.f90)
insert_sias_program(read_write_n00 ${SRC}/programs/read_write_n00.f90)
