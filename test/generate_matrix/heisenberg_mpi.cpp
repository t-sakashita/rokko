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
#include <rokko/eigen3.hpp>
#include <rokko/utility/heisenberg_hamiltonian.hpp>
#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>
#include <rokko/utility/math.hpp>
#include <rokko/collective.hpp>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  rokko::parallel_dense_ev solver;
  solver.initialize(argc, argv);
  rokko::grid g(MPI_COMM_WORLD);

  std::size_t L = 8;
  std::size_t num_bonds = L - 1;
  std::vector<std::pair<int, int>> lattice;
  for (std::size_t i=0; i<L-1; ++i) lattice.push_back(std::make_pair(i, i+1));

  int myrank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  constexpr int root = 0;

  if (myrank == root) {
    std::cout << "L=" << L << " num_bonds=" << num_bonds << std::endl;
    for (std::size_t i=0; i < num_bonds; ++i) {
      std::cout << lattice[i].first << " " << lattice[i].second << std::endl;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  const int p = rokko::find_power_of_two(nprocs);
  if (nprocs != (1 << p)) {
    throw std::invalid_argument("This program can be run only with 2^n MPI processes");
  }
  int N = 1 << (L-p);

  // creating column vectors which forms a heisenberg hamiltonian.
  int N_seq = 1 << L;
  Eigen::VectorXd recv_buffer(N_seq);
  Eigen::VectorXd buffer(N);
  Eigen::VectorXd v_seq(N_seq), w_seq(N_seq);
  Eigen::VectorXd v(N), w(N);
  for (int i=0; i<N; ++i) {
    // sequential version
    v_seq.setZero();
    w_seq.setZero();
    v_seq(i) = 1;
    if (myrank == root) {
      rokko::heisenberg_hamiltonian::multiply(L, lattice, v_seq.data(), w_seq.data());
      std::cout << "sequential version:" << std::endl;
      std::cout << w_seq.transpose() << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // MPI version
    MPI_Scatter(v_seq.data(), N, MPI_DOUBLE, v.data(), N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    w.setZero();
    rokko::heisenberg_hamiltonian::multiply(MPI_COMM_WORLD, L, lattice, v.data(), w.data(), buffer.data());
    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == myrank) {
        std::cout << "myrank=" << myrank << std::endl;
        std::cout << w.transpose() << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(w.data(), N, MPI_DOUBLE, recv_buffer.data(), N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0) {
      std::cout << "seq = " << w_seq.transpose() << std::endl;
      std::cout << "recv = " << recv_buffer.transpose() << std::endl;
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
  if (myrank == root) {
    rokko::heisenberg_hamiltonian::fill_diagonal(L, lattice, w_seq.data());
    std::cout << "fill_diagonal sequential version:" << std::endl;
    std::cout << "seq = " << w_seq.transpose() << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);    
  // MPI version
  if (myrank == root) {  
    std::cout << "fill_diagonal MPI version:" << std::endl;
  }
  rokko::heisenberg_hamiltonian::fill_diagonal(MPI_COMM_WORLD, L, lattice, w.data());
  for (int proc=0; proc<nprocs; ++proc) {
    if (proc == myrank) {
      std::cout << "myrank=" << myrank << std::endl;
      std::cout << w.transpose() << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // test for generate function
  Eigen::MatrixXd lmat(N_seq, N_seq);
  rokko::heisenberg_hamiltonian::generate(L, lattice, lmat);

  rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(N_seq, g);
  rokko::distributed_matrix<double,rokko::matrix_col_major> mat(map);
  rokko::heisenberg_hamiltonian::generate(L, lattice, mat);
  Eigen::MatrixXd lmat_gather(N_seq, N_seq);
  rokko::gather(mat, lmat_gather, root);

  if (myrank == root) {
    std:: cout << "lmat:" << std::endl << lmat << std::endl;
    std:: cout << "lmat_gather:" << std::endl << lmat_gather << std::endl;
    if (lmat_gather == lmat) {
      std::cout << "OK: distributed_matrix by 'generate' equals to a eigen_matrix by 'generate'." << std::endl;
    } else {
      std::cout << "ERROR: distributed_matrix by 'generate' is differnet from a eigen_matrix by 'generate'."<< std::endl;
      exit(1);
    }
  }

  MPI_Finalize();
  return 0;
}
