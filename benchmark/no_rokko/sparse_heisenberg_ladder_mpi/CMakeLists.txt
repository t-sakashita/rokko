#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2014-2014 by Synge Todo <wistaria@comp-phys.org>
#                            Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

if(BUILD_ANASAZI)
  set(SAMPLE anasazi_crs anasazi_mfree slepc_crs slepc_mfree anasazi_krylov-schur_crs anasazi_krylov-schur_mfree anasazi_crs_verify anasazi_mfree_verify)
  foreach(name ${SAMPLE})
    set(target_name benchmark_no_rokko_sparse_heisenberg_ladder_mpi_${name})
    add_executable(${target_name} ${name}.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko)
  endforeach(name)
endif(BUILD_ANASAZI)

