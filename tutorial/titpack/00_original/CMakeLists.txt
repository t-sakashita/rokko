# Copyright (C) 2014 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(titpack NONE)
enable_language(Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

add_library(titpack common.f large.f medium.f small.f clock.f90)

set(SAMPLES sample-01 sample-02 sample-03 sample-04 sample-05 sample-06 sample-07 sample-08 sample-09 sample-10 sample-11 sample-11-14)
foreach(name ${SAMPLES})
  add_executable(${name} ${name}.f)
  target_link_libraries(${name} titpack)
endforeach(name ${SAMPLES})
