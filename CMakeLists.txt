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
option(BUILD_EIGEN_S "Build eigen_s library" ON)
option(BUILD_EIGEN_SX "Build eigen_sx library" ON)

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
include_directories(${PROJECT_SOURCE_DIR}/3rd-party/eigen3)
if (BOOST_INCLUDE_DIR)
  include_directories(BEFORE ${BOOST_INCLUDE_DIR})
endif()

set(CMAKE_EXE_LINKER_FLAGS ${MPI_CXX_LINK_FLAGS})

# find Elemental
find_package(Elemental)

# test environment
enable_testing()
include(CTest)

# directories
add_subdirectory(3rd-party)
add_subdirectory(rokko)
add_subdirectory(sample_scalapack)
add_subdirectory(test)
add_subdirectory(sample_eigenK)
add_subdirectory(sample_elemental)
add_subdirectory(benchmark)
