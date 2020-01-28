/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_DISTRIBUTED_MFREE_TO_CRS_HPP
#define ROKKO_DISTRIBUTED_MFREE_TO_CRS_HPP

#include <rokko/distributed_crs_matrix.hpp>
#include <rokko/distributed_mfree.hpp>

namespace rokko {

void distributed_mfree_to_crs(rokko::distributed_mfree const& op, rokko::distributed_crs_matrix& mat) {
  if (mat.get_solver_name() == "anasazi") {
    throw std::invalid_argument("distributed_mfree_to_crs() : output_matrix_market is not available for Anasazi yet.");
  }

  std::vector<double> x(op.get_num_local_rows()), y(op.get_num_local_rows());
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  for (int global_row=0; global_row<op.get_dim(); ++global_row) {
    x.assign(op.get_num_local_rows(), 0.);
    y.assign(op.get_num_local_rows(), 0.);
    if ((global_row >= op.get_local_offset()) && (global_row < (op.get_local_offset() + op.get_num_local_rows()))) {
      x[global_row - op.get_local_offset()] = 1.;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    op.multiply(x.data(), y.data());
    MPI_Barrier(MPI_COMM_WORLD);
    for (int local_col=0; local_col<op.get_num_local_rows(); ++local_col) {
      if (y[local_col] != 0) {
        int global_col = local_col + op.get_local_offset();
        mat.insert(global_row, 1, &global_col, &y[local_col]);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  mat.complete();
}

} // end namespace rokko

#endif // ROKKO_DISTRIBUTED_MFREE_TO_CRS_HPP
