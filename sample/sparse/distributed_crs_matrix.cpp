/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 by Synge Todo <wistaria@comp-phys.org>,
*                       Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <mpi.h>
#include <iostream>

#include <rokko/parallel_sparse_solver.hpp>
#include <rokko/distributed_crs_matrix.hpp>

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  rokko::parallel_sparse_solver solver("slepc");
  int dim = 4;
  rokko::distributed_crs_matrix mat(dim, dim, solver);

  int num_nonzero_cols[] = {2, 1, 2, 1};
  int nonzero_cols[] = {0, 1, 3, 0, 3, 2};
  double values[] = {7.1, 5.2, 6.4, 0.2, 4.3, 0.5};

  int current = 0;
  for (int row = 0; row < dim; ++row) {
    std::cout << "current=" << current << std::endl;
    mat.insert(row, num_nonzero_cols[row], &nonzero_cols[current], &values[current]);
    current += num_nonzero_cols[row]; 
  }
  mat.complete();
  mat.print();

  MPI_Finalize();
}
