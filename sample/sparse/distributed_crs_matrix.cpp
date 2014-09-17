/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <mpi.h>
#include <iostream>

//#include <rokko/grid_1d.hpp>
#include <rokko/parallel_sparse_solver.hpp>
#include <rokko/distributed_crs_matrix.hpp>

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  int dim = 4;
  rokko::parallel_sparse_solver solver("anasazi");
  rokko::distributed_crs_matrix mat(dim, dim, solver);

  int rows[] = {1, 1, 2, 3, 3, 4};

  int num_nonzero_cols[] = {2, 1, 2, 1};
  int nonzero_cols[] = {0, 1, 3, 0, 3, 2};
  double values[] = {1., 2., 4., 1., 4., 3.};

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
