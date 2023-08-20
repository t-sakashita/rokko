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
#include <tuple>

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/xyz_hamiltonian.hpp>
#include <rokko/utility/xyz_hamiltonian_mpi.hpp>
#include <rokko/collective.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

TEST(xyz_hamiltonian, serial_mpi) {
  rokko::parallel_dense_ev solver;
  solver.initialize(global_argc, global_argv);
  const rokko::grid g(MPI_COMM_WORLD);

  constexpr std::size_t L = 4;
  constexpr auto num_bonds = L - 1;
  std::vector<std::pair<int, int>> lattice;
  std::vector<std::tuple<double, double, double>> coupling;
  for (std::size_t i=0; i<L-1; ++i) {
    lattice.emplace_back(std::make_pair(i, i+1));
    coupling.emplace_back(std::make_tuple(1, 0.3, 0.2));
  }

  int myrank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  constexpr int root = 0;

  if (myrank == root) {
    std::cout << "L=" << L << " num_bonds=" << num_bonds << std::endl;
    for (std::size_t i=0; i<num_bonds; ++i) {
      std::cout << lattice[i].first << " " << lattice[i].second << " "
                << std::get<0>(coupling[i]) << " " << std::get<1>(coupling[i]) << " " << std::get<2>(coupling[i]) << std::endl;
    }
  }

  const auto p = rokko::find_power_of_two(nprocs);
  if (nprocs != (1 << p)) {    
    throw std::invalid_argument("This program can be run only with 2^n MPI processes");
  }
  const auto N = 1 << (L-p);

  const auto N_seq = 1 << L;
  const auto N_seq_proc = (myrank == root) ? N_seq : 0;
  Eigen::VectorXd buffer(N);
  Eigen::VectorXd v_seq(N_seq_proc), w_seq(N_seq_proc), w_gather(N_seq_proc);
  Eigen::VectorXd v(N), w(N);

  // testing multiply
  for (int i=0; i<N; ++i) {
    // MPI version
    v.setZero();
    if (myrank == (i / N))  v(i % N) = 1;
    w.setZero();
    rokko::xyz_hamiltonian::multiply(MPI_COMM_WORLD, L, lattice, coupling, v.data(), w.data(), buffer.data());
    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == myrank) {
        std::cout << "myrank=" << myrank << std::endl;
        std::cout << w.transpose() << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Gather(w.data(), N, MPI_DOUBLE, w_gather.data(), N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    if (myrank == root) {
      // sequential version
      v_seq.setZero();
      v_seq(i) = 1;
      w_seq.setZero();
      rokko::xyz_hamiltonian::multiply(L, lattice, coupling, v_seq, w_seq);
      std::cout << "multiply: i=" << i << std::endl;
      std::cout << "seq = " << w_seq.transpose() << std::endl;
      std::cout << "gather = " << w_gather.transpose() << std::endl;
      ASSERT_TRUE(w_seq == w_gather);
    }
  }

  // testing fill_diagonal
  // MPI version
  rokko::xyz_hamiltonian::fill_diagonal(MPI_COMM_WORLD, L, lattice, coupling, w.data());
  for (int proc=0; proc<nprocs; ++proc) {
    if (proc == myrank) {
      std::cout << "myrank=" << myrank << std::endl;
      std::cout << w.transpose() << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Gather(w.data(), N, MPI_DOUBLE, w_gather.data(), N, MPI_DOUBLE, root, MPI_COMM_WORLD);
  if (myrank == root) {
    rokko::xyz_hamiltonian::fill_diagonal(L, lattice, coupling, w_seq.data());
    std::cout << "fill_diagonal:" << std::endl;
    std::cout << "seq = " << w_seq.transpose() << std::endl;
    std::cout << "gather = " << w_gather.transpose() << std::endl;
    ASSERT_TRUE(w_seq == w_gather);
  }

  // testing generate
  const rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(N_seq, g);
  rokko::distributed_matrix<double,rokko::matrix_col_major> mat(map);
  rokko::xyz_hamiltonian::generate(L, lattice, coupling, mat);
  Eigen::MatrixXd lmat_gather(N_seq_proc, N_seq_proc);
  rokko::gather(mat, lmat_gather, root);

  if (myrank == root) {
    Eigen::MatrixXd lmat(N_seq, N_seq);
    rokko::xyz_hamiltonian::generate(L, lattice, coupling, lmat);
    std::cout << "generate:" << std::endl;
    std::cout << "lmat:" << std::endl << lmat << std::endl;
    std::cout << "lmat_gather:" << std::endl << lmat_gather << std::endl;
    ASSERT_TRUE(lmat_gather == lmat);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  global_argc = argc;
  global_argv = argv;
  const auto result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
