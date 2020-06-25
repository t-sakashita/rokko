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

#include <mpi.h>
#include <cmath>
#include <iostream>
#include <rokko/distributed_matrix.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/elpa.hpp>
#include <rokko/elpa/diagonalize_set_parameters.hpp>

int main(int argc, char** argv) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  if (elpa_init(20200417) != ELPA_OK) {
    throw std::invalid_argument("ERROR: elpa::initialize()");
  }
  rokko::grid grid(MPI_COMM_WORLD);

  int n = 8;
  if (argc > 1) n = std::atoi(argv[1]);

  if (grid.get_myrank() == 0) {
    std::cout << "n = " << n << std::endl;
    std::cout << "nprocs = " << grid.get_nprocs() << std::endl;
    std::cout << "nprow = " << grid.get_nprow() << std::endl;
    std::cout << "npcol = " << grid.get_npcol() << std::endl;
  }

  rokko::mapping_bc<rokko::matrix_col_major> map(n, 1);
  rokko::distributed_matrix<double> a(map), eigvecs(map);
  Eigen::VectorXd eigvals(n);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      a.set_global(i, j, std::min(i, j) + 1);

  int error;
  elpa_t handle = elpa_allocate(&error);

  rokko::parameters params;
  params.set("routine", "elpa2");
  rokko::elpa::set_parameters(a, params, handle);
  assert_elpa_ok(elpa_setup(handle));

  rokko::elpa::set_solver(params, handle);
  int info = rokko::elpa::diag(handle, a, eigvals, eigvecs);
  assert_elpa_ok( rokko::elpa::deallocate(handle) );

  if (grid.get_myrank() == 0)
    std::cout << "eigenvalues: " << eigvals.transpose() << std::endl;

  MPI_Finalize();
  return 0;
}
