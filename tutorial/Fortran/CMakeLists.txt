# Copyright (C) 2014 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(wakate-handson NONE)
enable_language(CXX Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

find_package(Rokko REQUIRED PATHS ${ROKKO_ROOT_DIR} $ENV{ROKKO_HOME}
  NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found Rokko: ${ROKKO_ROOT_DIR}")
include(${ROKKO_USE_FILE})
message(STATUS ${ROKKO_USE_FILE})

add_executable(LAPACK LAPACK.f)
target_link_libraries(LAPACK ${ROKKO_LAPACK_LIBRARIES})

add_executable(ScaLAPACK ScaLAPACK.f)
#target_link_libraries(ScaLAPACK ${ROKKO_LIBRARIES})
target_link_libraries(ScaLAPACK ${ROKKO_SCALAPACK_LIBRARIES})


add_executable(Rokko Rokko.f90)
target_link_libraries(Rokko ${ROKKO_LIBRARIES}
  ${ROKKO_LIBRARIES_Fortran})
