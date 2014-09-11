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

// C++ version of TITPACK Ver.2 by H. Nishimori

/************ Sample main program #11 *****************
* 1d Heisenberg antiferromagnet with 8 spins
* Eigenvalues and an eigenvector by diag
******************************************************/

#include <mpi.h>
#include <iostream>
#include <boost/timer.hpp>
#include <rokko/solver.hpp>
#include "titpack.hpp"
#include "options.hpp"

typedef rokko::parallel_dense_solver solver_type;
typedef rokko::distributed_matrix<rokko::matrix_col_major> matrix_type;

int main(int argc, char** argv) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid g(comm, rokko::grid_col_major);
  int myrank = g.get_myrank();
  int root = 0;

  std::cout.precision(10);
  options opt(argc, argv, 8, solver_type::default_solver(), myrank == 0);
  if (!opt.valid) MPI_Abort(MPI_COMM_WORLD, 34);
  boost::timer tm;
  MPI_Barrier(comm);
  double t1 = tm.elapsed();

  // lattice structure
  int n = opt.N;
  int ibond = n;
  std::vector<int> ipair;
  for (int i = 0; i < ibond; ++i) {
    ipair.push_back(i);
    ipair.push_back((i + 1) % n);
  }

  // Hamiltonian parameters
  std::vector<double> bondwt(ibond, -1);
  std::vector<double> zrtio(ibond, 1);

  // table of configurations and Hamiltonian operator
  subspace ss(n, 0);
  hamiltonian hop(ss, ipair, bondwt, zrtio);
  solver_type solver(opt.solver);
  solver.initialize(argc, argv);

  // Hamiltonian matrix
  matrix_type elemnt(hop.dimension(), hop.dimension(), g, solver);
  elm3(hop, elemnt);
  MPI_Barrier(comm);
  double t2 = tm.elapsed();
  
  rokko::localized_vector E(hop.dimension());
  matrix_type v(hop.dimension(), hop.dimension(), g, solver);
  solver.diagonalize(elemnt, E, v);
  double t3 = tm.elapsed();

  if (myrank == 0) {
    int ne = 4;
    std::cout << "[Eigenvalues]\n";
    for (int i = 0; i < ne; ++i) std::cout << '\t' << E[i];
    std::cout << std::endl;
  }

  // // Do not forget to call elm3 again before calling check3
  elm3(hop, elemnt);
  // check3(elemnt, v, 0);
  MPI_Barrier(comm);
  double t4 = tm.elapsed();
  
  std::vector<int> npair;
  npair.push_back(1);
  npair.push_back(2);
  // std::vector<double> sxx(1);
  // xcorr(ss, npair, v, 0, sxx);
  // std::cout << "sxx: " << sxx[0] << std::endl;
  std::vector<double> szz(1);
  zcorr(comm, ss, npair, v, 0, szz);
  if (myrank == 0) std::cout << "szz: " << szz[0] << std::endl;
  MPI_Barrier(comm);
  double t5 = tm.elapsed();

  if (myrank == 0) {
    std::cerr << "initialize      " << (t2-t1) << " sec\n"
              << "diagonalization " << (t3-t2) << " sec\n"
              << "check           " << (t4-t3) << " sec\n"
              << "correlation     " << (t5-t4) << " sec\n";
  }
  MPI_Finalize();
}
