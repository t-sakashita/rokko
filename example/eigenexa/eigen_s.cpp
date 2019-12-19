/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <mpi.h>
#include <cmath>
#include <iostream>
#include <rokko/distributed_matrix.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/eigenexa.hpp>

int main(int argc, char** argv) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  rokko::grid grid(MPI_COMM_WORLD);
  rokko::eigenexa::init(grid.get_comm(), grid.is_row_major() ? 'R' : 'C');

  int n = 8;
  if (argc > 1) n = std::atoi(argv[1]);

  if (grid.get_myrank() == 0) {
    std::cout << "n = " << n << std::endl;
    std::cout << "nprocs = " << grid.get_nprocs() << std::endl;
    std::cout << "nprow = " << grid.get_nprow() << std::endl;
    std::cout << "npcol = " << grid.get_npcol() << std::endl;
  }

  auto nxy = rokko::eigenexa::get_matdims(grid, n);
  rokko::mapping_bc<rokko::matrix_col_major> map(n, 1, {nxy.first, nxy.second}, grid);
  rokko::distributed_matrix<double> a(map);
  rokko::distributed_matrix<double> z(map);
  Eigen::VectorXd w(n);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      a.set_global(i, j, std::min(i, j) + 1);
  
  rokko::eigenexa::eigen_s(a, w, z);
  if (grid.get_myrank() == 0)
    std::cout << "eigenvalues: " << w.transpose() << std::endl;

  rokko::eigenexa::free();
  MPI_Finalize();
}
