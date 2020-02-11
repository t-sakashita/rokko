/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_OUTPUT_MATRIX_MARKET_HPP
#define ROKKO_OUTPUT_MATRIX_MARKET_HPP

#include <rokko/distributed_mfree.hpp>
#include <rokko/utility/mpi_vector.hpp>
#include <iostream>

namespace rokko {

void output_matrix_market(distributed_mfree const& op, std::ostream& os = std::cout) {
  rokko::mpi_comm comm(op.get_comm());
  const int nprocs = comm.get_nprocs();
  const int myrank = comm.get_myrank();
  const int dim = op.get_dim();
  rokko::mpi_vector mpi(dim, op.get_comm());
  rokko::mapping_1d map(dim, comm);
  const int start_row = map.start_row();
  const int end_row = map.end_row();

  Eigen::Vector<double> x(op.get_num_local_rows()), y(op.get_num_local_rows());
  Eigen::Vector<double> vec(dim);

  constexpr int root_proc = 0;
  std::vector<int> cols;
  std::vector<double> values;

  if (comm.get_myrank() == root_proc) {
    os << "%%MatrixMarket matrix coordinate real general" << std::endl;
    os << dim << " " << dim << " " << 100 << std::endl; // Fix: nnz = 100
  }
  comm.barrier();

  for (int row=0; row<dim; ++row) {
    x.setZero();
    y.setZero();
    if ((row >= start_row) && (row < end_row)) {
      x[row - start_row] = 1.;
    }
    op.multiply(x.data(), y.data());
    mpi.gather(y, vec, root_proc);
    if (myrank == root_proc) {
      for (int col=0; col<dim; ++col) {
        if (vec[col] != 0) {
          os << row + 1 << " " << col+ 1 << " " << vec[col] << std::endl;
        }
      }
    }
  }
}

} // end namespace rokko

#endif // ROKKO_OUTPUT_MATRIX_MARKET_HPP
