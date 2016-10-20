/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string library(rokko::parallel_sparse_ev::default_solver());
  if (argc >= 2) library = argv[1];
  
  int dim = 8;
  int num_nonzero_cols[] = {2, 1, 1, 2, 3, 2, 2, 2};
  int nonzero_cols[] = {0, 4, 3, 5, 1, 7, 0, 5, 6, 2, 4, 4, 7, 3, 6};
  double values[] = {7.1, 2.8, 6.4, 0.5, 6.4, 3.5, 2.8, 0.2, 1.4, 0.5, 0.2, 1.4, 4.3, 3.5, 4.3};

  if (rank == 0) std::cout << "[solver = " << library << "]" << std::endl;
  rokko::parallel_sparse_ev solver(library);
  rokko::distributed_crs_matrix mat(dim, dim, solver);
  int current = 0;
  for (int row = 0; row < dim; ++row) {
    mat.insert(row, num_nonzero_cols[row], &nonzero_cols[current], &values[current]);
    current += num_nonzero_cols[row]; 
  }
  mat.complete();
  mat.print();
  
  solver.finalize();
  MPI_Finalize();
}
