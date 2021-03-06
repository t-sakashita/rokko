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

#include <rokko/solver.hpp>
#include <rokko/utility/mpi_vector.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

void run_test(MPI_Comm comm, int dim) {
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  rokko::mpi_vector mpi(dim);

  Eigen::VectorXd lvec(dim);
  lvec.setRandom();
#ifndef NDEBUG
  if (rank == 0) std::cout << lvec.transpose() << std::endl;
#endif

  for (int r = 0; r < size; ++r) {
    rokko::distributed_vector<double> vec(dim, mpi.get_displacement(), mpi.get_displacement() + mpi.get_count());
    mpi.scatter(lvec, vec, r);
    for (int i = 0; i < dim; ++i) {
      if (vec.has_global_index(i))
        ASSERT_EQ(vec.get_global(i), lvec(i));
    }
  }
}

TEST(mpi_communication, scatter_vector) {
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  int dim = 100;
  if (global_argc > 1) {
    dim = std::stoi(global_argv[1]);
  }
  if (rank == 0) std::cout << "dimension = " << dim << std::endl;

  run_test(comm, dim);
}

int main(int argc, char** argv) {
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  global_argc = argc;
  global_argv = argv;
  result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
