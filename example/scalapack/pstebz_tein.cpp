/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
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
#include <rokko/scalapack.hpp>
#include <rokko/utility/laplacian_matrix.hpp>

int main(int argc, char** argv) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  rokko::grid grid(MPI_COMM_WORLD);

  int n = 20;
  if (argc > 1) n = std::atoi(argv[1]);

  // generate matrix
  Eigen::VectorXd d(n), e(n-1);
  rokko::laplacian_matrix::generate(d, e);

  if (grid.get_myrank() == 0) {
    std::cout << "n = " << n << std::endl;
    std::cout << "nprocs = " << grid.get_nprocs() << std::endl;
    std::cout << "nprow = " << grid.get_nprow() << std::endl;
    std::cout << "npcol = " << grid.get_npcol() << std::endl;
    if ((grid.get_nprocs() != grid.get_nprow() * grid.get_npcol()) ||
        (n % grid.get_nprow() != 0) || (n % grid.get_npcol() != 0)) {
      std::cerr << "incompatible matrix size and number of processes\n";
      MPI_Abort(MPI_COMM_WORLD, 127);
    }
  }

  const int ictxt = grid.get_blacs_context();
  Eigen::VectorXd w(n);
  constexpr double abstol = 0.;
  int m, nsplit;
  Eigen::VectorXi iblock(n), isplit(n);

  int info = rokko::scalapack::pstebz(ictxt, 'I', 'E',
                                      0., 0., 0, 0,
                                      abstol, d, e, m, nsplit,
                                      w, iblock, isplit);

  rokko::mapping_bc<rokko::matrix_col_major> map(n, 1);
  rokko::distributed_matrix<double> z(map);
  double orfac = 0;
  Eigen::VectorXi ifail(map.get_m_global());
  Eigen::VectorXi iclustr(2 * map.get_nprow() * map.get_npcol());
  Eigen::VectorXd gap(map.get_nprow() * map.get_npcol());

  info = rokko::scalapack::pstein(d, e, m,
                                  w, iblock, isplit, orfac,
                                  z,
                                  ifail, iclustr, gap);

  if (info) {
    std::cerr << "Error in rokko::scalapack::psyevd\n";
    MPI_Abort(MPI_COMM_WORLD, info);
  }
  if (grid.get_myrank() == 0)
    std::cout << "eigenvalues: " << w.transpose() << std::endl;

  MPI_Finalize();
}
