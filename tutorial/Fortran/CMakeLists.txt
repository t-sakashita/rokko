#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(wakate-handson NONE)

find_package(Rokko REQUIRED PATHS ${ROKKO_ROOT_DIR} $ENV{ROKKO_HOME}
  NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found Rokko: ${ROKKO_ROOT_DIR}")
message(STATUS "Found path to use-rokko.cmake: ${ROKKO_USE_FILE}")
include(${ROKKO_USE_FILE})

add_executable(LAPACK LAPACK.f)
target_link_libraries(LAPACK ${ROKKO_LAPACK_LIBRARIES})

add_executable(ScaLAPACK ScaLAPACK.f)
#target_link_libraries(ScaLAPACK ${ROKKO_LIBRARIES})
target_link_libraries(ScaLAPACK ${ROKKO_SCALAPACK_LIBRARIES})


add_executable(Rokko Rokko.f90)
target_link_libraries(Rokko ${ROKKO_LIBRARIES})
