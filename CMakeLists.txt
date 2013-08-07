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

# enable C and C++ compilers
enable_language(C CXX Fortran)


# options
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_LAPACKE "Build lapacke library" ON)
option(BUILD_EIGEN_S "Build eigen_s library" ON)
option(BUILD_EIGEN_SX "Build eigen_sx library" ON)

# RPATH setting
if(APPLE)
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else(APPLE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif(APPLE)

# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  # Almost always OpenMP flags are same both for C and for Fortran.
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
endif(OPENMP_FOUND)

# MPI library
find_package(MPI)

set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(${PROJECT_SOURCE_DIR})
if (BUILD_LAPACKE)
  include_directories(${PROJECT_BINARY_DIR}/3rd-party/lapacke/include)
  include_directories(${PROJECT_SOURCE_DIR}/3rd-party/lapacke/include)
endif (BUILD_LAPACKE)
include_directories(${PROJECT_SOURCE_DIR}/3rd-party/eigen3)
if (BOOST_INCLUDE_DIR)
  include_directories(BEFORE ${BOOST_INCLUDE_DIR})
endif()

set(CMAKE_EXE_LINKER_FLAGS ${MPI_CXX_LINK_FLAGS})

# LAPACK/ScaLAPACK fallback
if ( NOT DEFINED SCALAPACK_LIB )
  # Find LAPACK libraries
  find_package(BLAS)
  if ( BLAS_FOUND )
    find_package(LAPACK)
    if ( LAPACK_FOUND )
      set(SCALAPACK_LIB ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      set(SCALAPACK_LIB ${SCALAPACK_LIB} -lscalapack-openmpi -lblacs-openmpi -lblacsCinit-openmpi)
    endif ( LAPACK_FOUND )
  endif ( BLAS_FOUND )
endif ( NOT DEFINED SCALAPACK_LIB )

if (EIGEN_EXA_LIB)
   set(BUILD_EIGEN_SX OFF)
   set(BUILD_EIGEN_EXA ON)
endif (EIGEN_EXA_LIB)

if (ELPA_LIB)
   set(BUILD_ELPA ON)
endif (ELPA_LIB)

# find Elemental
if (ELEMENTAL_DIR)
   find_package(Elemental)
endif (ELEMENTAL_DIR)

# find PETSc/SLEPc
# Normally PETSc is built with MPI, if not, use CC=mpicc, etc, so we need to specify CXX as an option.
if (USE_SLEPC)
   set(PETSC_DIR "/home/sakashita/rokko_lib")
   set(SLEPC_DIR "/home/sakashita/rokko_lib")
   find_package(PETSc2 COMPONENTS CXX)
   message(STATUS "Display_main PETSC_LIBRARIES=${PETSC_LIBRARIES}.")
   # Added by Sakashita
   set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDES})
   message(STATUS "Display_main PETSC_INCLUDE_DIRS=${PETSC_INCLUDE_DIRS}.")
   find_package(SLEPc2)
   message(STATUS "Display_main CMAKE_REQUIRED_INCLUDES=${CMAKE_REQUIRED_INCLUDES}.")
   message(STATUS "Display_main CMAKE_REQUIRED_LIBRARIES=${CMAKE_REQUIRED_LIBRARIES}.")
endif (USE_SLEPC)

# test environment
enable_testing()
include(CTest)

# directories
add_subdirectory(3rd-party)
add_subdirectory(rokko)
add_subdirectory(sample)
add_subdirectory(sample_scalapack)
add_subdirectory(test)
add_subdirectory(sample_eigenK)
add_subdirectory(sample_elemental)
add_subdirectory(benchmark)
add_subdirectory(c-fortran_sample)
add_subdirectory(sample_slepc)

#set(TRILINOS_PREFIX "/home/sakashita/build/install_build")
if( DEFINED TRILINOS_PREFIX )
  add_subdirectory(sample_anasazi)
endif( DEFINED TRILINOS_PREFIX )
