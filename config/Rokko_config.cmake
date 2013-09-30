#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
#                            Synge Todo <wistaria@comp-phys.org>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# options
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_LAPACKE "Build lapacke library" ON)

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

find_package(BLAS)
find_package(LAPACK)
find_package(SCALAPACK)
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
find_package(PETSc2 COMPONENTS CXX)
find_package(SLEPc2)
set(BUILD_PETSC ${PETSC_FOUND})
set(BUILD_SLEPC ${SLEPC_FOUND})
