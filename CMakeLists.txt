#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
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

# RPATH setting
if(APPLE)
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else(APPLE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif(APPLE)

# enable C and C++ compilers
enable_language(C CXX Fortran)

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
find_package(MPI)
set(CMAKE_EXE_LINKER_FLAGS ${MPI_CXX_LINK_FLAGS})

set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR}/3rd-party/lapacke/include)
include_directories(${PROJECT_SOURCE_DIR}/3rd-party/lapacke/include)
include_directories(${PROJECT_SOURCE_DIR}/3rd-party/eigen3)
if (BOOST_INCLUDE_DIR)
  include_directories(BEFORE ${BOOST_INCLUDE_DIR})
endif()

# find_package(BLAS)
find_package(Lapack)
if(NOT LAPACK_FOUND)
  message(FATAL_ERROR "LAPACK not found")
endif(NOT LAPACK_FOUND)

# find Scalapack
find_package(Scalapack)
set(BUILD_SCALAPACK ${SCALAPACK_FOUND})

# find EigenExa
find_package(EigenExa)
set(BUILD_EIGENEXA ${EIGENEXA_FOUND})

# find Elemental
find_package(Elemental)
set(BUILD_ELEMENTAL ${ELEMENTAL_FOUND})

# find Elpa
find_package(Elpa)
set(BUILD_ELPA ${ELPA_FOUND})

# find PETSc/SLEPc
find_package(PETSc COMPONENTS CXX)
find_package(SLEPc)
set(BUILD_PETSC ${PETSC_FOUND})
set(BUILD_SLEPC ${SLEPC_FOUND})


configure_file(${PROJECT_SOURCE_DIR}/config/RokkoConfig.cmake.in ${PROJECT_BINARY_DIR}/config/UseRokko.cmake @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/config/add_rokko_test.cmake ${PROJECT_BINARY_DIR}/config/add_barista_test.cmake @COPYONLY)
#configure_file(${PROJECT_SOURCE_DIR}/config/add_rokko_test.cmake ${PROJECT_BINARY_DIR}/config/add_barista_test.cmake @COPYONLY)
#configure_file(${PROJECT_SOURCE_DIR}/config/add_rokko_test.cmake ${PROJECT_BINARY_DIR}/config/add_barista_test.cmake @COPYONLY)

install(FILES ${PROJECT_SOURCE_DIR}/config/add_rokko_test.cmake  DESTINATION config/)
install(FILES ${PROJECT_SOURCE_DIR}/config/run_test.cmake  DESTINATION config/)
install(FILES ${PROJECT_SOURCE_DIR}/config/run_test_mpi.cmake  DESTINATION config/)


# test environment
enable_testing()
include(CTest)

# directories
add_subdirectory(3rd-party)
add_subdirectory(rokko)
add_subdirectory(sample)
add_subdirectory(test)
add_subdirectory(benchmark)
add_subdirectory(c-fortran_sample)
add_subdirectory(sample_slepc)

#set(TRILINOS_PREFIX "/home/sakashita/build/install_build")
if( DEFINED TRILINOS_PREFIX )
  add_subdirectory(sample_anasazi)
endif( DEFINED TRILINOS_PREFIX )
