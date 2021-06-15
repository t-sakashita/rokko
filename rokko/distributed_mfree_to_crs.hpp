/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <rokko/distributed_crs_matrix.hpp>
#include <rokko/distributed_mfree.hpp>
#include <rokko/utility/mpi_vector.hpp>

namespace rokko {

void distributed_mfree_to_crs(rokko::distributed_mfree const& op, rokko::distributed_crs_matrix& mat) {
  rokko::mpi_comm mpi_comm(op.get_comm());
  const int nprocs = mpi_comm.get_nprocs();
  const int dim = op.get_dim();
  rokko::mpi_vector mpi(dim, op.get_comm());
  rokko::mapping_1d map(dim, mpi_comm);
  const int start_row = map.start_row();
  const int end_row = map.end_row();

  Eigen::Vector<double> x(op.get_num_local_rows()), y(op.get_num_local_rows());
  Eigen::Vector<double> vec(dim);

  for (int row=0; row<dim; ++row) {
    x.setZero();
    y.setZero();
    if ((row >= start_row) && (row < end_row)) {
      x[row - start_row] = 1.;
    }
    op.multiply(x.data(), y.data());
    const int proc = (nprocs * row) / dim;
    mpi.gather(y, vec, proc);
    if ((row >= start_row) && (row < end_row)) {
      for (int col=0; col<dim; ++col) {
        if (vec[col] != 0) {
          mat.insert(row, 1, &col, &vec[col]);
        }
      }
    }
  }
  mat.complete();
}

} // end namespace rokko
