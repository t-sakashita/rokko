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

include(Rokko_config)

configure_file(${PROJECT_SOURCE_DIR}/config/RokkoConfig.cmake.in ${PROJECT_BINARY_DIR}/config/UseRokko.cmake @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/config/add_rokko_test.cmake ${PROJECT_BINARY_DIR}/config/add_barista_test.cmake @COPYONLY)
#configure_file(${PROJECT_SOURCE_DIR}/config/add_rokko_test.cmake ${PROJECT_BINARY_DIR}/config/add_barista_test.cmake @COPYONLY)
#configure_file(${PROJECT_SOURCE_DIR}/config/add_rokko_test.cmake ${PROJECT_BINARY_DIR}/config/add_barista_test.cmake @COPYONLY)

install(FILES ${PROJECT_SOURCE_DIR}/config/add_rokko_test.cmake  DESTINATION config/)
install(FILES ${PROJECT_SOURCE_DIR}/config/run_test.cmake  DESTINATION config/)
install(FILES ${PROJECT_SOURCE_DIR}/config/run_test_mpi.cmake  DESTINATION config/)


# test environment
enable_testing()
include(CTest)

# directories
add_subdirectory(3rd-party)
add_subdirectory(rokko)
add_subdirectory(sample)
add_subdirectory(test)
add_subdirectory(benchmark)
add_subdirectory(c-fortran_sample)
add_subdirectory(sample_slepc)

#set(TRILINOS_PREFIX "/home/sakashita/build/install_build")
if( DEFINED TRILINOS_PREFIX )
  add_subdirectory(sample_anasazi)
endif( DEFINED TRILINOS_PREFIX )
