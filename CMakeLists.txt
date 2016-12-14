#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

cmake_minimum_required(VERSION 2.8)

project(rokko NONE)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)
include(add_rokko_test)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_C_BINDING "Build C binding libraries" ON)
option(BUILD_Fortran_BINDING "Build Fortran binding libraries" ON)
option(BUILD_Python_BINDING "Build Python binding libraries" OFF)
option(DISABLE_OPENMP "Disable OpenMP Parallelization" OFF)
option(DISABLE_MPI "Disable MPI Parallelization" OFF)
option(DISABLE_SCALAPACK "Disable ScaLAPACK Library" OFF)
option(DISABLE_EIGENEXA "Disable EigenExa Library" OFF)
option(DISABLE_ELEMENTAL "Disable Elemental Library" OFF)
option(DISABLE_ELPA "Disable ELPA Library" OFF)
option(DISABLE_PETSC "Disable PETSc Library" OFF)
option(DISABLE_SLEPC "Disable SLEPc Library" OFF)
option(DISABLE_ANASAZI "Disable Anasazi Library" OFF)

option(DISABLE_PDSYEVR "DISABLE pdsyevr of ScaLAPACK" OFF)

option(BUILD_TESTS "Build test suites" ON)
option(BUILD_EXAMPLES "Build examples" ON)
option(BUILD_TUTORIALS "Build tutorials" ON)
option(BUILD_BENCHMARKS "Build benchmark tests" OFF)
option(BUILD_DEVELOPERS "Build codes for Rokko's developers" OFF)

if(NOT BUILD_C_BINDING)
  set(BUILD_Fortran_BINDING OFF)
endif(NOT BUILD_C_BINDING)

if(DISABLE_MPI)
  set(DISABLE_SCALAPACK ON)
  set(DISABLE_EIGENEXA ON)
  set(DISABLE_ELEMENTAL ON)
  set(DISABLE_ELPA ON)
  set(DISABLE_PETSC ON)
  set(DISABLE_SLEPC ON)
  set(DISABLE_ANASAZI ON)
endif(DISABLE_MPI)

if(DISABLE_SCALAPACK)
  set(DISABLE_EIGENEXA ON)
endif(DISABLE_SCALAPACK)
if(DISABLE_PETSC)
  set(DISABLE_SLEPC ON)
endif(DISABLE_PETSC)

# RPATH setting
if(APPLE AND (CMAKE_VERSION VERSION_LESS 2.8.12))
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else(APPLE AND (CMAKE_VERSION VERSION_LESS 2.8.12))
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  set(CMAKE_MACOSX_RPATH 1)
endif(APPLE AND (CMAKE_VERSION VERSION_LESS 2.8.12))

# C++, C, and Fortran compilers
enable_language(CXX C Fortran)
set(ROKKO_HAVE_C_BINDING ${BUILD_C_BINDING})
if(BUILD_Fortran_BINDING)
#  enable_language(Fortran)
  # ancient gfortran does not support iso_c_bindings
  try_compile(CMake_Fortran_SUPPORTS_ISO_C_BINDING
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/config/check_iso_c_binding.f90
    OUTPUT_VARIABLE LOG)
  if(NOT CMake_Fortran_SUPPORTS_ISO_C_BINDING)
    set(BUILD_Fortran_BINDING OFF)
    MESSAGE(STATUS "Do not build Fortran binding since the fortran compiler does not support iso_c_binding.")
  endif(NOT CMake_Fortran_SUPPORTS_ISO_C_BINDING)
endif(BUILD_Fortran_BINDING)
set(ROKKO_HAVE_Fortran_BINDING ${BUILD_Fortran_BINDING})

# OpenMP
if(NOT DISABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    # Almost always OpenMP flags are same both for C and for Fortran.
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
    set(ALPS_USE_MKL_PARALLEL TRUE)
  else(OPENMP_FOUND)
    set(DISABLE_OPENMP ON)
  endif(OPENMP_FOUND)
endif(NOT DISABLE_OPENMP)
if(NOT DISABLE_OPENMP)
  set(ROKKO_HAVE_OPENMP True)
else(NOT DISABLE_OPENMP)
  set(ROKKO_HAVE_OPENMP False)
endif(NOT DISABLE_OPENMP)

# MPI library
if(NOT DISABLE_MPI)
  find_package(MPI)
  if(MPI_FOUND)
    set(CMAKE_EXE_LINKER_FLAGS ${MPI_CXX_LINK_FLAGS})
  else(MPI_FOUND)
    set(DISABLE_MPI ON)
  endif(MPI_FOUND)
endif(NOT DISABLE_MPI)
if(NOT DISABLE_MPI)
  set(ROKKO_HAVE_MPI True)
else(NOT DISABLE_MPI)
  set(ROKKO_HAVE_MPI False)
endif(NOT DISABLE_MPI)

# Linker for Fortran Program
option(USE_CXX_LINKER "Use C++ linker for Fortran program" OFF)
message(STATUS "Use C++ linker for Fortran program: " ${USE_CXX_LINKER})
if(USE_CXX_LINKER)
  set(ROKKO_USE_CXX_LINKER 1)
endif(USE_CXX_LINKER)

# standard include directories
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR}/3rd-party/lapacke/include)
include_directories(${PROJECT_SOURCE_DIR}/3rd-party/lapacke/include)
include_directories(${PROJECT_SOURCE_DIR}/3rd-party/eigen3)

# Python
if(BUILD_Python_BINDING)
  find_package(PythonInterp REQUIRED)
  find_package(PythonLibs REQUIRED)
  include_directories(${PYTHON_INCLUDE_DIRS})
endif(BUILD_Python_BINDING)

# Boost Library
if(BUILD_Python_BINDING)
  find_package(Boost COMPONENTS python system)
else(BUILD_Python_BINDING)
  find_package(Boost COMPONENTS system)
endif(BUILD_Python_BINDING)
include_directories(BEFORE ${Boost_INCLUDE_DIR})

if((EXISTS $ENV{MKLROOT}) AND (NOT BLA_VENDOR))
  set(BLA_VENDOR Intel10_64lp)
endif((EXISTS $ENV{MKLROOT}) AND (NOT BLA_VENDOR))
find_package(LAPACK)
if(NOT LAPACK_FOUND)
  message(FATAL_ERROR "LAPACK not found")
endif(NOT LAPACK_FOUND)

#
# parallel dense solvers
#

# find ScaLAPACK
if(DISABLE_SCALAPACK)
  set(BUILD_SCALAPACK FALSE)
else(DISABLE_SCALAPACK)
  find_package(ScaLAPACK)
  set(BUILD_SCALAPACK ${SCALAPACK_FOUND})
endif(DISABLE_SCALAPACK)
set(ROKKO_HAVE_SCALAPACK ${BUILD_SCALAPACK})
string(REGEX REPLACE ";" " " SCAPALACK_LIBS "${SCAPALACK_LIBRARIES}")

if(DISABLE_PDSYEVR)
  set(ROKKO_HAVE_PDSYEVR OFF)
else(DISABLE_PDSYEVR)
  set(ROKKO_HAVE_PDSYEVR ON)
endif(DISABLE_PDSYEVR)

# find EigenExa
if(DISABLE_EIGENEXA)
  set(BUILD_EIGENEXA FALSE)
else(DISABLE_EIGENEXA)
  find_package(EigenExa)
  set(BUILD_EIGENEXA ${EIGENEXA_FOUND})
  include_directories(${EIGENEXA_INCLUDE_DIR})
endif(DISABLE_EIGENEXA)
set(ROKKO_HAVE_EIGENEXA ${BUILD_EIGENEXA})
string(REGEX REPLACE ";" " " EIGENEXA_LIBS "${EIGENEXA_LIBRARIES}")

# find Elemental
if(DISABLE_ELEMENTAL)
  set(BUILD_ELEMENTAL FALSE)
else(DISABLE_ELEMENTAL)
  find_package(Elemental)
  set(BUILD_ELEMENTAL ${ELEMENTAL_FOUND})
  include_directories(${ELEMENTAL_INCLUDE_DIR})
endif(DISABLE_ELEMENTAL)
set(ROKKO_HAVE_ELEMENTAL ${BUILD_ELEMENTAL})
string(REGEX REPLACE ";" " " ELEMENTAL_LIBS "${ELEMENTAL_LIBRARIES}")
if(ELEMENTAL_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11")
  message(STATUS "Compiler option: ${CMAKE_CXX_FLAGS}")
endif(ELEMENTAL_FOUND)

# find ELPA
if(DISABLE_ELPA)
  set(BUILD_ELPA FALSE)
else(DISABLE_ELPA)
  find_package(ELPA)
  set(BUILD_ELPA ${ELPA_FOUND})
endif(DISABLE_ELPA)
set(ROKKO_HAVE_ELPA ${BUILD_ELPA})
string(REGEX REPLACE ";" " " ELPA_LIBS "${ELPA_LIBRARIES}")

if(BUILD_SCALAPACK OR BUILD_EIGENEXA OR BUILD_ELEMENTAL OR BUILD_ELPA)
  set(BUILD_PARALLEL_DENSE_SOLVER ON)
else(BUILD_SCALAPACK OR BUILD_EIGENEXA OR BUILD_ELEMENTAL OR BUILD_ELPA)
  set(BUILD_PARALLEL_DENSE_SOLVER OFF)
endif(BUILD_SCALAPACK OR BUILD_EIGENEXA OR BUILD_ELEMENTAL OR BUILD_ELPA)
set(ROKKO_HAVE_PARALLEL_DENSE_SOLVER ${BUILD_PARALLEL_DENSE_SOLVER})

#
# parallel sparse solvers
#

# find PETSc/SLEPc
if(DISABLE_PETSC)
  set(BUILD_PETSC On)
else(DISABLE_PETSC)
  find_package(PETSc COMPONENTS CXX)
  set(BUILD_PETSC ${PETSC_FOUND})
  include_directories(${PETSC_INCLUDE_DIR})
endif(DISABLE_PETSC)
if(DISABLE_SLEPC)
  set(BUILD_SLEPC OFF)
else(DISABLE_SLEPC)
  find_package(SLEPc)
  set(BUILD_SLEPC ${SLEPC_FOUND})
  include_directories(${SLEPC_INCLUDE_DIR})
endif(DISABLE_SLEPC)
set(ROKKO_HAVE_PETSC ${BUILD_PETSC})
set(ROKKO_HAVE_SLEPC ${BUILD_SLEPC})
string(REGEX REPLACE ";" " " PETSC_LIBS "${PETSC_LIBRARIES}")
string(REGEX REPLACE ";" " " SLEPC_LIBS "${SLEPC_LIBRARIES}")

# find Anasazi
if(DISABLE_ANASAZI)
  set(BUILD_ANASAZI FALSE)
else(DISABLE_ANASAZI)
  find_package(Anasazi)
  set(BUILD_ANASAZI ${ANASAZI_FOUND})
  include_directories(${ANASAZI_INCLUDE_DIR})
endif(DISABLE_ANASAZI)
set(ROKKO_HAVE_ANASAZI ${BUILD_ANASAZI})
string(REGEX REPLACE ";" " " ANASAZI_LIBS "${ANASAZI_LIBRARIES}")

if(BUILD_SLEPC OR BUILD_ANASAZI)
  set(BUILD_PARALLEL_SPARSE_SOLVER ON)
else(BUILD_SLEPC OR BUILD_ANASAZI)
  set(BUILD_PARALLEL_SPARSE_SOLVER OFF)
endif(BUILD_SLEPC OR BUILD_ANASAZI)
set(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER ${BUILD_PARALLEL_SPARSE_SOLVER})

set(ROKKO_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")

configure_file(${PROJECT_SOURCE_DIR}/rokko/config.h.in ${PROJECT_BINARY_DIR}/rokko/config.h)
install(FILES ${PROJECT_BINARY_DIR}/rokko/config.h COMPONENT headers DESTINATION include/rokko/)
configure_file(${PROJECT_SOURCE_DIR}/config/RokkoConfig.cmake.in ${PROJECT_BINARY_DIR}/config/RokkoConfig.cmake @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/config/include.mk.in ${PROJECT_BINARY_DIR}/config/include.mk @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/config/RokkoConfig.cmake ${PROJECT_BINARY_DIR}/config/include.mk ${PROJECT_SOURCE_DIR}/config/UseRokko.cmake ${PROJECT_SOURCE_DIR}/config/add_rokko_test.cmake ${PROJECT_SOURCE_DIR}/config/run_test.cmake ${PROJECT_SOURCE_DIR}/config/run_test_mpi.cmake DESTINATION share/rokko/)

# enable C++11
if(BUILD_ELPA)
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif(BUILD_ELPA)

# test environment
enable_testing()
include(CTest)

# directories
add_subdirectory(3rd-party)
add_subdirectory(rokko)
add_subdirectory(tool)

if(BUILD_TESTS)
  add_subdirectory(test)
endif(BUILD_TESTS)

if(BUILD_EXAMPLES)
  add_subdirectory(example)
  if(BUILD_DEVELOPERS)
    add_subdirectory(developer)
  endif(BUILD_DEVELOPERS)
endif(BUILD_EXAMPLES)

if(BUILD_BENCHMARKS)
  add_subdirectory(benchmark)
endif(BUILD_BENCHMARKS)

if(BUILD_TUTORIALS)
  add_subdirectory(tutorial)
endif(BUILD_TUTORIALS)
