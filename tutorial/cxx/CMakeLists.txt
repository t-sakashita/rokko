#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(myproject_cxx NONE)

find_package(Rokko REQUIRED PATHS ${ROKKO_ROOT_DIR} $ENV{ROKKO_HOME}
  NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found Rokko: ${ROKKO_ROOT_DIR}")
include(${ROKKO_USE_FILE})  # informing compilers, include and library paths, CMAKE_BUILD_TYPE and the supported language bindings

display_rokko_components()

add_executable(minij_mpi minij_mpi.cpp)
target_link_libraries(minij_mpi ${ROKKO_LIBRARIES})
