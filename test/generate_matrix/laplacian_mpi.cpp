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

#include <rokko/rokko.hpp>
#include <rokko/utility/laplacian_matrix.hpp>
#include <rokko/utility/laplacian_mfree.hpp>
#include <rokko/utility/mpi_vector.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

TEST(laplacian_mfree, serial_mpi) {
  int myrank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  constexpr int root = 0;

  constexpr int N_seq = 20;
  const int N_seq_proc = (myrank == root) ? N_seq : 0;
  Eigen::VectorXd v_seq(N_seq_proc), w_seq(N_seq_proc), recv_buffer(N_seq_proc);

  rokko::laplacian_mfree op(N_seq);
  const int N = op.get_num_local_rows();
  Eigen::VectorXd v(N), w(N);

  rokko::mpi_vector mpi(N_seq);

  for (int i=0; i<N_seq; ++i) {
    if (myrank == root) {
      v_seq.setZero();
      v_seq(i) = 1;
      w_seq.setZero();
      rokko::laplacian_matrix::multiply(N_seq, v_seq.data(), w_seq.data());
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // MPI version
    mpi.scatter(v_seq, v);
    w.setZero();
    op.multiply(v.data(), w.data());
    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == myrank) {
        std::cout << "myrank=" << myrank << std::endl;
        std::cout << w.transpose() << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    mpi.gather(w, recv_buffer);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == root) {
      std::cout << "i=" << i << std::endl;
      std::cout << "w_seq=" << w_seq.transpose() << std::endl;
      std::cout << "recv=" << recv_buffer.transpose() << std::endl;
      ASSERT_TRUE(w_seq == recv_buffer);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

int main(int argc, char** argv) {
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  global_argc = argc;
  global_argv = argv;
  result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
