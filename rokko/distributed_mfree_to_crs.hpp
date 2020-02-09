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

#ifndef ROKKO_DISTRIBUTED_MFREE_TO_CRS_HPP
#define ROKKO_DISTRIBUTED_MFREE_TO_CRS_HPP

#include <rokko/distributed_crs_matrix.hpp>
#include <rokko/distributed_mfree.hpp>
#include <rokko/utility/mpi_vector.hpp>

namespace rokko {

void distributed_mfree_to_crs(rokko::distributed_mfree const& op, rokko::distributed_crs_matrix& mat) {
  MPI_Comm comm = op.get_comm();
  rokko::mpi_comm mpi_comm(comm);
  const int nprocs = mpi_comm.get_nprocs();
  const int myrank = mpi_comm.get_myrank();
  const int dim = op.get_dim();
  rokko::mpi_vector mpi(dim, comm);
  rokko::mapping_1d map(dim, mpi_comm);

  Eigen::Vector<double> x(op.get_num_local_rows()), y(op.get_num_local_rows());
  Eigen::Vector<double> vec(dim);

  for (int row=0; row<dim; ++row) {
    x.setZero();
    y.setZero();
    if ((row >= op.get_local_offset()) && (row < (op.get_local_offset() + op.get_num_local_rows()))) {
      x[row - op.get_local_offset()] = 1.;
    }
    MPI_Barrier(comm);
    op.multiply(x.data(), y.data());
    MPI_Barrier(comm);
    const int proc = (nprocs * row) / dim;
    mpi.gather(y, vec, proc);
    if ((row >= map.start_row()) && (row < map.end_row())) {
      for (int col=0; col<dim; ++col) {
        if (vec[col] != 0) {
          mat.insert(row, 1, &col, &vec[col]);
        }
      }
    }
    MPI_Barrier(comm);
  }
  MPI_Barrier(comm);
  mat.complete();
}

} // end namespace rokko

#endif // ROKKO_DISTRIBUTED_MFREE_TO_CRS_HPP
