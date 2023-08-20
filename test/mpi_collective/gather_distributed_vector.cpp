/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2019 Rokko Developers https://github.com/t-sakashita/rokko
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

template <typename T, int MAJOR=Eigen::ColMajor>
using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1, MAJOR>;

template <typename T>
void run_test(MPI_Comm comm, int dim) {
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  rokko::mpi_vector mpi(dim);

  rokko::distributed_vector<T> vec(dim, mpi.get_displacement(), mpi.get_displacement() + mpi.get_count());
  for (int i = 0; i < dim; ++i) {
    if (vec.has_global_index(i)) vec.set_global(i, static_cast<T>(i));
  }

  VectorX<T> lvec(dim);
  for (int r = 0; r < size; ++r) {
    mpi.gather(vec, lvec, r);
    if (rank == r) {
#ifndef NDEBUG
      std::cout << "rank=" << rank << std::endl
                << lvec.transpose() << std::endl;
      std::cout.flush();
#endif
      for (int i = 0; i < dim; ++i) {
        if (vec.has_global_index(i))
          ASSERT_EQ(vec.get_global(i), static_cast<T>(i));
      }
    }
#ifndef NDEBUG
    MPI_Barrier(comm);
#endif
  }

}

TEST(mpi_communication, gather_vector) {
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  const int dim = (global_argc > 1) ? std::stoi(global_argv[1]) : 100;
  if (rank == 0) std::cout << "dimension = " << dim << std::endl;

  run_test<float>(comm, dim);
  run_test<double>(comm, dim);
  run_test<long double>(comm, dim);
  run_test<int>(comm, dim);
  //run_test<long int>(comm, dim);
  run_test<char>(comm, dim);
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
