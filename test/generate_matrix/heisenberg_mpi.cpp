/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <iostream>
#include <vector>

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/utility/heisenberg_hamiltonian.hpp>
#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>
#include <rokko/collective.hpp>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  rokko::parallel_dense_ev solver;
  solver.initialize(argc, argv);
  rokko::grid g(MPI_COMM_WORLD);

  int L = 8;
  int num_bonds = L - 1;
  std::vector<std::pair<int, int>> lattice;
  for (int i=0; i<L-1; ++i) lattice.push_back(std::make_pair(i, i+1));

  int myrank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Status status;
  const int root = 0;
  int ierr;

  if (myrank == root) {
    std::cout << "L=" << L << " num_bonds=" << num_bonds << std::endl;
    for (int i=0; i < num_bonds; ++i) {
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

  if (nprocs != (1 << p)) {
    if ( myrank == 0 ) {
      std::cout << "This program can be run only for powers of 2" << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  int N = 1 << (L-p);

  // creating column vectors which forms a heisenberg hamiltonian.
  int N_seq = 1 << L;
  double* recv_buffer = new double[N_seq];
  double* send_buffer = new double[N];
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
    v.assign(N, -2.);
    MPI_Scatter(&v_seq[0], N, MPI_DOUBLE, &v[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&w[0], N, MPI_DOUBLE, recv_buffer, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0) {
      std::cout << "seq=";
      for (int j=0; j<N_seq; ++j) {
	std::cout << w_seq[j] << " ";
      }
      std::cout << std::endl;
      std::cout << "recv=";
      for (int j=0; j<N_seq; ++j) {
	std::cout << recv_buffer[j] << " ";
      }
      std::cout << std::endl;
      for (int j=0; j<N_seq; ++j) {
	if (w_seq[j] != recv_buffer[j]) {
	  std::cout << "j=" << j << "  w_seq[j]=" << w_seq[j] << "  recv[j]=" << recv_buffer[j] << std::endl;
	  exit(1);
	}
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

  }

  // test fill_diagonal of quantum heisenberg hamiltonian.
  // sequential version
  std::vector<double> w_seq;
  w_seq.assign(N_seq, 0);
  if (myrank == root) {
    rokko::heisenberg_hamiltonian::fill_diagonal(L, lattice, w_seq);
    std::cout << "fill_diagonal sequential version:" << std::endl;
    for (int j=0; j<N_seq; ++j) {
      std::cout << w_seq[j] << " ";
    }
    std::cout << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);    
  // MPI version
  std::vector<double> w;
  w.assign(N, 0);
  if (myrank == root) {  
    std::cout << "fill_diagonal MPI version:" << std::endl;
  }
  rokko::heisenberg_hamiltonian::fill_diagonal(MPI_COMM_WORLD, L, lattice, w);
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

  // test for generate function
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> lmat(N_seq, N_seq);
  rokko::heisenberg_hamiltonian::generate(L, lattice, lmat);

  rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(N_seq, g);
  rokko::distributed_matrix<double,rokko::matrix_col_major> mat(map);
  rokko::heisenberg_hamiltonian::generate(L, lattice, mat);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> lmat_gather(N_seq, N_seq);
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
