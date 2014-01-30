#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
#                            Synge Todo <wistaria@comp-phys.org>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

cmake_minimum_required(VERSION 2.8)

project(rokko NONE)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)
include(add_rokko_test)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_C_BINDING "Build C binding libraries" ON)
option(BUILD_Fortran_BINDING "Build Fortran binding libraries" ON)
option(DISABLE_SCALAPACK "Disable Scalapack Library" OFF)
option(DISABLE_EIGENEXA "Disable EigenExa Library" OFF)
option(DISABLE_ELEMENTAL "Disable Elemental Library" OFF)
option(DISABLE_ELPA "Disable Elpa Library" OFF)
option(DISABLE_PETSC "Disable PETSc Library" OFF)
option(DISABLE_SLEPC "Disable SLEPc Library" OFF)
option(DISABLE_ANASAZI "Disable Anasazi Library" OFF)

option(DISABLE_MPI "Disable MPI support" OFF)
option(DISABLE_CUDA "Disable CUDA support" ON)

if(DISABLE_SCALAPACK)
  set(DISABLE_EIGENEXA ON)
endif(DISABLE_SCALAPACK)
if(DISABLE_PETSC)
  set(DISABLE_SLEPC ON)
endif(DISABLE_PETSC)

# MPI disable settings
if(DISABLE_MPI)
  set(DISABLE_SCALAPACK ON)
  set(DISABLE_EIGENEXA ON)
  set(DISABLE_ELEMENTAL ON)
  set(DISABLE_ELPA ON)
  set(DISABLE_SLEPC ON)
  set(BUILD_C_BINDING OFF)
  set(BUILD_Fortran_BINDING OFF)
endif(DISABLE_MPI)

# RPATH setting
if(APPLE)
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else(APPLE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif(APPLE)

# C++, C, and Fortran compilers
enable_language(CXX)
enable_language(C)
if(BUILD_C_BINDING)
  enable_language(C)
endif(BUILD_C_BINDING)
if(BUILD_Fortran_BINDING)
  enable_language(Fortran)
  # ancient gfortran does not support iso_c_bindings
  try_compile(CMake_Fortran_SUPPORTS_ISO_C_BINDING
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/config/check_iso_c_binding.f90
    OUTPUT_VARIABLE LOG)
  #MESSAGE("LOG: ${LOG}")
  if(NOT CMake_Fortran_SUPPORTS_ISO_C_BINDING)
    set(BUILD_Fortran_BINDING OFF)
    MESSAGE("Do not build Fortran binding since the fortran compiler does not support iso_c_binding.")
  endif(NOT CMake_Fortran_SUPPORTS_ISO_C_BINDING)
endif(BUILD_Fortran_BINDING)

# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  # Almost always OpenMP flags are same both for C and for Fortran.
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
  set(ALPS_USE_MKL_PARALLEL TRUE)
endif(OPENMP_FOUND)

# MPI library
if(DISABLE_MPI)
  set(BUILD_MPI OFF)
else(DISABLE_MPI)
  find_package(MPI)
  if(MPI_FOUND)
    set(CMAKE_EXE_LINKER_FLAGS ${MPI_CXX_LINK_FLAGS})
    set(BUILD_MPI ON)
    option(HAVE_MPI "Compile with MPI support" ON)
  else(MPI_FOUND)
    set(BUILD_MPI OFF)
  endif(MPI_FOUND)
endif(DISABLE_MPI)

# Boost Library
find_package(BoostSrc)

set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR}/3rd-party/lapacke/include)
include_directories(${PROJECT_SOURCE_DIR}/3rd-party/lapacke/include)
include_directories(${PROJECT_SOURCE_DIR}/3rd-party/eigen3)
include_directories(BEFORE ${Boost_INCLUDE_DIR})

# find_package(BLAS)
find_package(Lapack)
if(NOT LAPACK_FOUND)
  message(FATAL_ERROR "LAPACK not found")
endif(NOT LAPACK_FOUND)
if(MAC_VECLIB)
  set(LAPACK_LIBRARIES "-framework vecLib")
endif(MAC_VECLIB) 

# find Scalapack
if(DISABLE_SCALAPACK)
  set(BUILD_SCALAPACK FALSE)
else(DISABLE_SCALAPACK)
  find_package(Scalapack)
  set(BUILD_SCALAPACK ${SCALAPACK_FOUND})
endif(DISABLE_SCALAPACK)

# find EigenExa
if(DISABLE_EIGENEXA)
  set(BUILD_EIGENEXA FALSE)
else(DISABLE_EIGENEXA)
  find_package(EigenExa)
  set(BUILD_EIGENEXA ${EIGENEXA_FOUND})
endif(DISABLE_EIGENEXA)

# find Elemental
if(DISABLE_ELEMENTAL)
  set(BUILD_ELEMENTAL FALSE)
else(DISABLE_ELEMENTAL)
  find_package(Elemental)
  set(BUILD_ELEMENTAL ${ELEMENTAL_FOUND})
endif(DISABLE_ELEMENTAL)

# find Elpa
if(DISABLE_ELPA)
  set(BUILD_ELPA FALSE)
else(DISABLE_ELPA)
  find_package(Elpa)
  set(BUILD_ELPA ${ELPA_FOUND})
endif(DISABLE_ELPA)

# find PETSc/SLEPc
if(DISABLE_PETSC)
  set(BUILD_PETSC FALSE)
else(DISABLE_PETSC)
  find_package(PETSc COMPONENTS CXX)
  set(BUILD_PETSC ${PETSC_FOUND})
endif(DISABLE_PETSC)
if(DISABLE_SLEPC)
  set(BUILD_SLEPC FALSE)
else(DISABLE_SLEPC)
  find_package(SLEPc)
  set(BUILD_SLEPC ${SLEPC_FOUND})
endif(DISABLE_SLEPC)

# find Anasazi
if(DISABLE_ANASAZI)
  set(BUILD_ANASAZI FALSE)
else(DISABLE_ANASAZI)
  find_package(Anasazi)
  set(BUILD_ANASAZI ${ANASAZI_FOUND})
endif(DISABLE_ANASAZI)

# find CUDA
if(DISABLE_CUDA)
  set(BUILD_CUDA FALSE)
else(DISABLE_CUDA)
  enable_language(Fortran)
  find_package(CUDA)
  option(HAVE_CUBLAS "cublas support" 1)
  if(CUDA_FOUND)
    find_package(MAGMA)
    set(BUILD_MAGMA ${MAGMA_FOUND})
    option(HAVE_MAGMA "MAGMA support" ${MAGMA_FOUND})
    include_directories(${MAGMA_INCLUDE_DIRS})
    include_directories(${CUDA_TOOLKIT_INCLUDE})
  endif(CUDA_FOUND)
endif(DISABLE_CUDA)

set(ROKKO_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/include")
set(ROKKO_LIBRARIES "${CMAKE_INSTALL_PREFIX}/lib/librokko.so")
set(ROKKO_FORTRAN_LIBRARIES "${CMAKE_INSTALL_PREFIX}/lib/librokko_fortran.so")

configure_file(${PROJECT_SOURCE_DIR}/config/RokkoConfig.cmake.in ${PROJECT_BINARY_DIR}/config/UseRokko.cmake @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/config/add_rokko_test.cmake ${PROJECT_BINARY_DIR}/config/add_barista_test.cmake @COPYONLY)

install(FILES ${PROJECT_SOURCE_DIR}/config/add_rokko_test.cmake  DESTINATION config/)
install(FILES ${PROJECT_SOURCE_DIR}/config/run_test.cmake  DESTINATION config/)
install(FILES ${PROJECT_SOURCE_DIR}/config/run_test_mpi.cmake  DESTINATION config/)
install(FILES ${PROJECT_BINARY_DIR}/config/UseRokko.cmake  DESTINATION config/)

# config file
configure_file(
  "${PROJECT_SOURCE_DIR}/rokko/config.hpp.in"
  "${PROJECT_BINARY_DIR}/rokko/config.hpp"
)

# Add the binary tree to the search path for include files
# so that we will find config.hpp
include_directories("${PROJECT_BINARY_DIR}")

# test environment
enable_testing()
include(CTest)

# directories
add_subdirectory(3rd-party)
add_subdirectory(rokko)
if (BUILD_Fortran_BINDING AND BUILD_MPI)
  add_subdirectory(rokko_c-fortran)
endif (BUILD_Fortran_BINDING AND BUILD_MPI)

if(BUILD_MPI)
add_subdirectory(sample_dense)
endif(BUILD_MPI)
add_subdirectory(test)
if(BUILD_MPI)
add_subdirectory(benchmark_dense)
endif(BUILD_MPI)
add_subdirectory(sample_slepc)
add_subdirectory(sample_anasazi)
add_subdirectory(benchmark_anasazi)


if (BUILD_Fortran_BINDING AND BUILD_MPI)
  add_subdirectory(sample_c)
  add_subdirectory(sample_fortran)
endif (BUILD_Fortran_BINDING AND BUILD_MPI)

