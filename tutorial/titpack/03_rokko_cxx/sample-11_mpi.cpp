/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// C++ version of TITPACK Ver.2 by H. Nishimori

/************ Sample main program #11 *****************
* 1d Heisenberg antiferromagnet with 8 spins
* Eigenvalues and an eigenvector by diag
******************************************************/

#include <mpi.h>
#include <iostream>
#include <chrono>
#include <rokko/rokko.hpp>
#include "titpack.hpp"
#include "options.hpp"

using solver_type = rokko::parallel_dense_ev;
using matrix_major = rokko::matrix_col_major;
using matrix_type = rokko::distributed_matrix<double, matrix_major>;

int main(int argc, char** argv) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  rokko::grid g;

  std::cout.precision(10);
  options opt(argc, argv, 8, solver_type::default_solver(), g.get_myrank() == 0);
  if (!opt.valid) MPI_Abort(MPI_COMM_WORLD, 1);
  MPI_Barrier(g.get_comm());
  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

  // lattice structure
  int n = opt.N;
  int ibond = n;
  std::vector<int> ipair;
  for (int i = 0; i < ibond; ++i) {
    ipair.emplace_back(i);
    ipair.emplace_back((i + 1) % n);
  }

  // Hamiltonian parameters
  std::vector<double> bondwt(ibond, -1);
  std::vector<double> zrtio(ibond, 1);

  // table of configurations and Hamiltonian operator
  subspace ss(n, 0);
  hamiltonian hop(ss, ipair, bondwt, zrtio);
  solver_type solver(opt.solver);
  rokko::mapping_bc<matrix_major> map = solver.default_mapping(hop.dimension(), g);

  solver.initialize(argc, argv);

  // Hamiltonian matrix
  matrix_type elemnt(map);
  elm3(hop, elemnt);
  MPI_Barrier(g.get_comm());
  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
  
  Eigen::VectorXd E(hop.dimension());
  matrix_type v(map);
  solver.diagonalize(elemnt, E, v);
  std::chrono::system_clock::time_point t3 = std::chrono::system_clock::now();

  if (g.get_myrank() == 0) {
    int ne = 4;
    std::cout << "[Eigenvalues]\n";
    for (int i = 0; i < ne; ++i) std::cout << '\t' << E[i];
    std::cout << std::endl << std::flush;
  }
  MPI_Barrier(g.get_comm());

  // Do not forget to call elm3 again before calling check3
  elm3(hop, elemnt);
  matrix_type w(map);
  check3_mpi(elemnt, v, 0, w);
  std::cout << std::flush;
  MPI_Barrier(g.get_comm());
  std::chrono::system_clock::time_point t4 = std::chrono::system_clock::now();
  
  std::vector<int> npair;
  npair.emplace_back(1);
  npair.emplace_back(2);
  std::vector<double> sxx(1);
  matrix_type sxmat(map);
  xcorr3_mpi(ss, npair, v, 0, sxx, sxmat, w);
  if (v.has_global_indices({0, 0})) std::cout << "sxx: " << sxx[0] << std::endl;
  std::cout << std::flush;
  MPI_Barrier(g.get_comm());
  std::vector<double> szz(1);
  zcorr_mpi(ss, npair, v, 0, szz);
  if (g.get_myrank() == 0) std::cout << "szz: " << szz[0] << std::endl;
  std::cout << std::flush;
  MPI_Barrier(g.get_comm());
  std::chrono::system_clock::time_point t5 = std::chrono::system_clock::now();

  if (g.get_myrank() == 0) {
    std::cerr << "initialize      " << 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count() << " sec\n"
              << "diagonalization " << 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(t3-t2).count() << " sec\n"
              << "check           " << 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(t4-t3).count() << " sec\n"
              << "correlation     " << 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(t5-t4).count() << " sec\n";
  }
  MPI_Finalize();
}
