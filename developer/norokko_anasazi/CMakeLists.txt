#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2014-2014 by Synge Todo <wistaria@comp-phys.org>
#                            Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

#include(/Users/sakashitatatsuya/rokko/trilinos-11.12.1-1/Debug/lib/cmake/Anasazi/AnasaziConfig.cmake)

if(BUILD_ANASAZI)
  set(SAMPLE heisenberg_crs_mpi heisenberg_mfree_mpi laplacian_mfree_mpi TraceMinDavidsonUserOpEx laplacian_mfree_parameters_mpi)
  foreach(name ${SAMPLE})
    set(target_name norokko_anasazi_${name})
    add_executable(${target_name} ${name}.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko)
  endforeach(name)
endif(BUILD_ANASAZI)

if(BUILD_ANASAZI)
  set(TESTS xyz_crs_mpi xyz_matfree_mpi mvtest)
  foreach(name ${TESTS})
    set(target_name norokko_anasazi_${name})
    add_executable(${target_name} ${name}.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko)
  endforeach(name)
endif(BUILD_ANASAZI)
