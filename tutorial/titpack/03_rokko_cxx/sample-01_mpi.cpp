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

/************* Sample main program #1 *****************
* 1d Heisenberg antiferromagnet with 16 spins
* Eigenvalues and an eigenvector / lnc1, lncv1
* Precision check and correlation functions
******************************************************/

#include <mpi.h>
#include <iostream>
#include <chrono>
#include <rokko/rokko.hpp>
#include "titpack.hpp"
#include "options.hpp"

using solver_type = rokko::parallel_sparse_ev;

int main(int argc, char** argv) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  rokko::mpi_comm comm;
  if (comm.get_nprocs() > 1) {
    std::cerr << "currently works only for Np = 1\n";
  }

  std::cout.precision(10);
  options opt(argc, argv, 16, solver_type::default_solver(), comm.get_myrank() == 0);
  if (!opt.valid) MPI_Abort(MPI_COMM_WORLD, 1);
  MPI_Barrier(comm.get_comm());
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
  solver.initialize(argc, argv);
  MPI_Barrier(comm.get_comm());
  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
  
  // Eigenvalues
  rokko::parameters params;
  params.set("Block Size", 5);
  params.set("Maximum Iterations", 500);
  params.set("Convergence Tolerance", 1.0e-8);
  params.set("num_eigenvalues", 10);
  solver.diagonalize(hop, params);
  std::chrono::system_clock::time_point t3 = std::chrono::system_clock::now();
  
  if (comm.get_myrank() == 0) {
    std::cout << "[Number of converged eigenpairs]\n\t" << solver.num_conv() << std::endl;
    // std::cout << "[Iteration number]\n\t" << itr << std::endl;
    std::cout << "[Eigenvalues]\n";
    for (int i = 0; i < solver.num_conv(); ++i) std::cout << '\t' << solver.eigenvalue(i);
    std::cout << std::endl;
  }

  // Ground-state eigenvector
  rokko::distributed_vector<double> eigvec;
  solver.eigenvector(0, eigvec);
  if (comm.get_myrank() == 0) std::cout << "[Eigenvector components (selected)]";
  std::cout << std::flush;
  MPI_Barrier(comm.get_comm());
  int count = 0;
  for (int i = 12; i < ss.dimension(); i += ss.dimension()/20, ++count) {
    if (eigvec.is_gindex(i)) {
      if (count % 4 == 0) std::cout << std::endl;
      std::cout << '\t' << eigvec.get_global(i);
    }
    std::cout << std::flush;
    MPI_Barrier(comm.get_comm());
  }
  if (comm.get_myrank() == 0) std::cout << std::endl;
  std::cout << std::flush;
  MPI_Barrier(comm.get_comm());

  // Precision check and correlation functions
  // double Hexpec = check2(mat, x, 0, v, 0);

  // std::vector<int> npair;
  // npair.emplace_back(1);
  // npair.emplace_back(2);
  // std::vector<double> sxx(1), szz(1);
  // xcorr(ss, npair, x, 0, sxx);
  // zcorr(ss, npair, x, 0, szz);
  // std::cout << "[Nearest neighbor correlation functions]\n\t" 
  //           << "sxx : " << sxx[0]
  //           << ", szz : " << szz[0] << std::endl;

  if (comm.get_myrank() == 0) {
    std::cerr << "initialize      " << 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count() << " sec\n"
              << "diagonalization " << 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(t3-t2).count() << " sec\n";
    // << "check           " << 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(t4-t3).count() << " sec\n"
    // << "correlation     " << 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(t5-t4).count() << " sec\n";
  }

  MPI_Finalize();
}
