#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2014-2014 by Synge Todo <wistaria@comp-phys.org>
#                            Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

if(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)
  include_directories(${ANASAZI_INCLUDE_DIR})
  set(SAMPLE heisenberg_crs_mpi heisenberg_mfree_mpi distributed_crs_matrix laplacian_crs_mpi laplacian_mfree_mpi laplacian_crs_parameter_mpi laplacian_mfree_parameter_mpi)
  foreach(name ${SAMPLE})
    set(target_name sample_sparse_${name})
    add_executable(${target_name} ${name}.cpp)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko)
  endforeach(name)
endif(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)

