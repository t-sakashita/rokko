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
  
  std::vector<std::string> solvers;
  if (argc >= 2) {
    solvers.push_back(argv[1]);
  } else {
    solvers = rokko::parallel_sparse_ev::solvers();
  }

  int dim = 4;
  int num_nonzero_cols[] = {2, 1, 2, 1};
  int nonzero_cols[] = {0, 1, 3, 0, 3, 2};
  double values[] = {7.1, 5.2, 6.4, 0.2, 4.3, 0.5};
  BOOST_FOREACH(std::string const& name, solvers) {
    if (rank == 0) std::cout << "[solver = " << name << "]" << std::endl;
    rokko::parallel_sparse_ev solver(name);
    rokko::distributed_crs_matrix mat(dim, dim, solver);
    int current = 0;
    for (int row = 0; row < dim; ++row) {
      mat.insert(row, num_nonzero_cols[row], &nonzero_cols[current], &values[current]);
      current += num_nonzero_cols[row]; 
    }
    mat.complete();
    mat.print();
  }
  
  MPI_Finalize();
}
