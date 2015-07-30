#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2014-2014 by Synge Todo <wistaria@comp-phys.org>
#                            Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

if(BUILD_SLEPC)
  set(TESTS frank_matfree_1proc heisenberg_crs_mpi heisenberg_matfree_1proc heisenberg_matfree_mpi xyz_crs_mpi xyz_matfree_mpi frank_matfree heisenberg_matfree_mpi_allvec)

  foreach(name ${TESTS})
    set(target_name norokko_slepc_diagonalize_${name})
    add_executable(${target_name} ${name}.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko)
  endforeach(name)
endif(BUILD_SLEPC)
