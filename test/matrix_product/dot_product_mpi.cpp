/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <rokko/collective.hpp>

#include <gtest/gtest.h>

constexpr double eps = 1e-11;

int global_argc;
char** global_argv;

TEST(dot_product_mpi, dot_product_mpi) {
  constexpr int root_proc = 0;
  const MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  const int dim = (global_argc > 1) ? std::stoi(global_argv[1]) : 100;

  if (rank == root_proc) std::cout << "dimension = " << dim << std::endl;
  const rokko::grid g(comm);
  const rokko::mapping_bc<rokko::matrix_col_major> map({dim, 1}, {1, 1}, g);
  rokko::distributed_matrix<double, rokko::matrix_col_major> vecX(map);
  rokko::distributed_matrix<double, rokko::matrix_col_major> vecY(map);
  Eigen::VectorXd locX(dim), locY(dim);
  locX.setRandom();
  locY.setRandom();
  rokko::scatter(locX, vecX, root_proc);
  rokko::scatter(locY, vecY, root_proc);

  // local calculation
  const double product_local = locX.dot(locY);

  // global calculation
  const double product_global = rokko::dot_product(vecX, false, 0, vecY, false, 0);

  for (int i = 0; i < dim; ++i) {
    if (vecX.has_global_indices({i, 0})) {
      std::cerr << i << ' ' << rank << ' ' << product_local << ' ' << product_global << std::endl;
      std::cerr.flush();
      EXPECT_NEAR(product_local, product_global, eps);
    }
    MPI_Barrier(MPI_COMM_WORLD);
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
