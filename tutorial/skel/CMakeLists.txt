#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(myproject NONE)

find_package(Rokko REQUIRED HINTS ${ROKKO_ROOT_DIR} $ENV{ROKKO_ROOT})
message(STATUS "Found Rokko: ${ROKKO_ROOT_DIR}")
include(${ROKKO_USE_FILE})
display_rokko_components()

enable_testing()
include(CTest)

add_executable(syev syev.cpp)
target_link_libraries(syev ${ROKKO_LIBRARIES})
add_test(syev syev)

add_executable(dsyev dsyev.c)
target_link_libraries(dsyev ${ROKKO_LIBRARIES})
add_test(dsyev dsyev)
