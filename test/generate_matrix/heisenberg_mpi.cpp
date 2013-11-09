/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/tuple/tuple.hpp>

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/heisenberg_hamiltonian.hpp>
#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>

#include <rokko/collective.hpp>

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  std::string solver_name("scalapack");
  rokko::solver solver(solver_name);
  solver.initialize(argc, argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid g(comm, rokko::grid_col_major);
  const int root = 0;

  /*
  if (argc <= 1) {
    std::cerr << "error: " << argv[0] << " path_to_heisenberg.ip" << std::endl;
    exit(1);
  }

  std::ifstream ifs(argv[1]);
  if (!ifs) {
    std::cout << "can't open file" << std::endl;
    exit(2);
  }
  */

  int L, num_bonds;
  std::vector<std::pair<int, int> > lattice;

  L = 4;
  num_bonds = L - 1;
  /*
  ifs >> L >> num_bonds;
  for (int i=0; i<num_bonds; ++i) {
    int j, k;
    ifs >> j >> k;
    lattice.push_back(std::make_pair(j, k));
  }
  
  for (int i=0; i<num_bonds; ++i) {
    double jx, jy, jz;
    ifs >> jx >> jy >> jz;
    coupling.push_back(boost::make_tuple(jx, jy, jz));
  }
  */
  for (int i=0; i<L-1; ++i) {
    lattice.push_back(std::make_pair(i, i+1));
  }

  int myrank, nprocs;
  MPI_Status status;
  int ierr;

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (myrank == root) {
    std::cout << "L=" << L << " num_bonds=" << num_bonds << std::endl;
    for (int i=0; i<num_bonds; ++i) {
      std::cout << lattice[i].first << " " << lattice[i].second << std::endl;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  int n = nprocs;
  int p = -1;
  do {
    n /= 2;
    ++p;
  } while (n > 0);

  //std::cerr << "nproc=" << nproc << " p=" << p << std::endl;
  if (nprocs != (1 << p)) {    
    if ( myrank == 0 ) {
      std::cout << "This program can be run only for powers of 2" << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  int N = 1 << (L-p);

  // creating column vectors which forms a heisenberg hamiltonian.
  int N_seq = 1 << L;
  std::vector<double> buffer(N);
  for (int i=0; i<N; ++i) {
    // sequential version
    std::vector<double> v_seq, w_seq;
    v_seq.assign(N_seq, 0);
    v_seq[i] = 1;
    w_seq.assign(N_seq, 0);
    if (myrank == root) {
      rokko::heisenberg_hamiltonian::multiply(L, lattice, v_seq, w_seq);
      std::cout << "sequential version:" << std::endl;
      for (int j=0; j<N_seq; ++j) {
        std::cout << w_seq[j] << " ";
      }
      std::cout << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // MPI version
    std::vector<double> v, w;
    v.assign(N, 0);
    if (myrank == (i / N))
      v[i % N] = 1;
    w.assign(N, 0);
    rokko::heisenberg_hamiltonian::multiply(MPI_COMM_WORLD, L, lattice, v, w, buffer);
    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == myrank) {
        std::cout << "myrank=" << myrank << std::endl;
        for (int j=0; j<N; ++j) {
          std::cout << w[j] << " ";
        }
        std::cout << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  // test for generate function
  rokko::localized_matrix<rokko::matrix_col_major> lmat(N_seq, N_seq);
  rokko::heisenberg_hamiltonian::generate(L, lattice, lmat);
  
  rokko::distributed_matrix<rokko::matrix_col_major> mat(N_seq, N_seq, g, solver);
  rokko::heisenberg_hamiltonian::generate(L, lattice, mat);
  rokko::localized_matrix<rokko::matrix_col_major> lmat_gather(N_seq, N_seq);
  rokko::gather(mat, lmat_gather, root);

  if (myrank == root) {
    std:: cout << "lmat:" << std::endl << lmat << std::endl;
    std:: cout << "lmat_gather:" << std::endl << lmat_gather << std::endl;
    if (lmat_gather == lmat) {
      std::cout << "OK: distributed_matrix by 'generate' equals to a localized_matrix by 'generate'." << std::endl;
    } else {
      std::cout << "ERROR: distributed_matrix by 'generate' is differnet from a localized_matrix by 'generate'."<< std::endl;
      exit(1);
    }
  }

  MPI_Finalize();
  return 0;
}


