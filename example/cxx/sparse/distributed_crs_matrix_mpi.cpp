/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
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

  const std::string library = (argc >= 2) ? argv[1] : rokko::parallel_sparse_ev::default_solver();

  constexpr int dim = 4;
  const int num_nonzero_cols[] = {2, 1, 2, 1};
  const int cols[] = {0, 1, 3, 0, 3, 2};
  const double values[] = {7.1, 5.2, 6.4, 0.2, 4.3, 0.5};

  if (rank == 0) std::cout << "[solver = " << library << "]" << std::endl;
  rokko::parallel_sparse_ev solver(library);
  rokko::distributed_crs_matrix mat(solver.default_mapping(dim, rokko::mpi_comm{MPI_COMM_WORLD}), 2);
  int current = 0;
  for (int row = 0; row < dim; ++row) {
    mat.insert(row, num_nonzero_cols[row], &cols[current], &values[current]);
    current += num_nonzero_cols[row];
  }
  mat.complete();
  mat.print();

  solver.finalize();
  MPI_Finalize();
}
