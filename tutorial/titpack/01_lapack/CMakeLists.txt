#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(titpack NONE)

find_package(Rokko REQUIRED PATHS ${ROKKO_ROOT_DIR} $ENV{ROKKO_HOME}
  NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found Rokko: ${ROKKO_ROOT_DIR}")
message(STATUS "Found path to use-rokko.cmake: ${ROKKO_USE_FILE}")
include(${ROKKO_USE_FILE})

add_library(titpack common.f large.f medium.f small.f clock.f90)

set(SAMPLES sample-01 sample-02 sample-03 sample-04 sample-05 sample-06 sample-07 sample-08 sample-09 sample-10 sample-11 sample-11-14)
foreach(name ${SAMPLES})
  add_executable(${name} ${name}.f)
  target_link_libraries(${name} titpack ${ROKKO_LIBRARIES})
  SET_TARGET_PROPERTIES(${name} PROPERTIES LINKER_LANGUAGE CXX)
  add_rokko_test(${name})
endforeach(name ${SAMPLES})
