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
#include <boost/lexical_cast.hpp>
#include <random>

#include <gtest/gtest.h>

constexpr double eps = 1e-11;

int global_argc;
char** global_argv;

TEST(product_v_mpi, product_v_mpi) {
  constexpr int root_proc = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  int dim = 100;
  if (global_argc > 1) {
    dim = boost::lexical_cast<int>(global_argv[1]);
  }

  std::mt19937 engine(123lu);
  std::uniform_real_distribution<> dist(-1.0, 1.0);

  if (rank == root_proc) std::cout << "dimension = " << dim << std::endl;
  rokko::parallel_dense_ev solver(rokko::parallel_dense_ev::default_solver());
  rokko::grid g(comm);
  rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(dim, g);
  rokko::distributed_matrix<double, rokko::matrix_col_major> matA(map);
  rokko::mapping_bc<rokko::matrix_col_major> mapvec({dim, 1}, {1, 1}, g);
  rokko::distributed_matrix<double, rokko::matrix_col_major> vecX(mapvec);
  rokko::distributed_matrix<double, rokko::matrix_col_major> vecY(mapvec);
  Eigen::MatrixXd locA(dim, dim);
  Eigen::VectorXd locX(dim), locY(dim);
  for (int j = 0; j < dim; ++j) for (int i = 0; i < dim; ++i) locA(i, j) = dist(engine);
  for (int i = 0; i < dim; ++i) locX(i) = dist(engine);
  for (int i = 0; i < dim; ++i) locY(i) = dist(engine);
  rokko::scatter(locA, matA, root_proc);
  rokko::scatter(locX, vecX, root_proc);
  rokko::scatter(locY, vecY, root_proc);

  // local calculation
  locY += locA * locX;

  // global calculation
  rokko::product_v(1.0, matA, false, vecX, false, 0, 1, vecY, false, 0);

  for (int i = 0; i < dim; ++i) {
    if (vecY.is_gindex({i, 0}))
      EXPECT_NEAR(vecY.get_global(i,0), locY(i), eps);
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
