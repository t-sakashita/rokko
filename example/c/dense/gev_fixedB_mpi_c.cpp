/*****************************************************************************
 *
 * Rokko: Integrated Interface for libraries of eigenvalue decomposition
 *
 * Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
 *                            Synge Todo <wistaria@comp-phys.org>,
 *                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
 *    
 * Distributed under the Boost Software License, Version 1.0. (See accompanying
 * file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *****************************************************************************/

#include <mpi.h>
#include <rokko/rokko.h>
#include <rokko/rokko_dense.h>

#include <stdio.h>
#include <stdlib.h>

#include "gev_fixedB_mpi.h"
#include "gev_fixedB_mpi.hpp"

typedef rokko::matrix_col_major matrix_major;

void diagonalize_fixedB_c(struct rokko_parallel_dense_solver solver_in, struct rokko_distributed_matrix A_in, struct rokko_distributed_matrix B_in,
			  struct rokko_localized_vector eigval_in, struct rokko_distributed_matrix eigvec_in, double tol) {
  rokko::parallel_dense_solver& solver = *(static_cast<rokko::parallel_dense_solver*>(solver_in.ptr));
  rokko::distributed_matrix<double, matrix_major>& A = *(static_cast<rokko::distributed_matrix<double,matrix_major>*>(A_in.ptr));
  rokko::distributed_matrix<double, matrix_major>& B = *(static_cast<rokko::distributed_matrix<double,matrix_major>*>(B_in.ptr));
  rokko::distributed_matrix<double, matrix_major>& eigvec = *(static_cast<rokko::distributed_matrix<double,matrix_major>*>(eigvec_in.ptr));
  rokko::localized_vector<double>& eigval = *(static_cast<rokko::localized_vector<double>*>(eigval_in.ptr));

  rokko::mapping_bc<matrix_major> map = A.get_mapping();
  rokko::distributed_matrix<double, matrix_major> tmp(map), Binvroot(map), mat(map);

  diagonalize_fixedB(solver, A, B, eigval, eigvec);
}

void set_A_B_c(struct rokko_localized_matrix locA_in, struct rokko_localized_matrix locB_in) {
  rokko::localized_matrix<double, matrix_major>& locA = *(static_cast<rokko::localized_matrix<double,matrix_major>*>(locA_in.ptr));
  rokko::localized_matrix<double, matrix_major>& locB = *(static_cast<rokko::localized_matrix<double,matrix_major>*>(locB_in.ptr));
  set_A_B(locA, locB);
}
