# Copyright (C) 2014 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(titpack NONE)
enable_language(CXX Fortran)

find_package(ROKKO REQUIRED PATHS ${ROKKO_ROOT_DIR} $ENV{ROKKO_HOME}
  NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ROKKO: ${ROKKO_ROOT_DIR}")
include(${ROKKO_USE_FILE})
message(STATUS ${ROKKO_USE_FILE})

add_library(titpack common.cpp large.cpp small.cpp) # medium.cpp

# set(SAMPLES sample-01 sample-02 sample-03 sample-04 sample-05 sample-06 sample-07 sample-08 sample-09 sample-10 sample-11)
set(SAMPLES sample-01 sample-02 sample-03 sample-11 checklists)
foreach(name ${SAMPLES})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} titpack ${ROKKO_LIBRARIES})
  add_rokko_test(${name})
endforeach(name ${SAMPLES})
