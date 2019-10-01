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

int global_argc;
char** global_argv;

TEST(dot_product_mpi, dot_product_mpi) {
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  int dim = 100;
  if (global_argc > 1) {
    dim = boost::lexical_cast<int>(global_argv[1]);
  }

  std::mt19937 engine(123lu);
  std::uniform_real_distribution<> dist(-1.0, 1.0);

  if (rank == 0) std::cout << "dimension = " << dim << std::endl;
  rokko::grid g(comm);
  rokko::mapping_bc<rokko::matrix_col_major> map(dim, 1, g, 1, 1);
  rokko::distributed_matrix<double, rokko::matrix_col_major> vecX(map);
  rokko::distributed_matrix<double, rokko::matrix_col_major> vecY(map);
  rokko::localized_matrix<double, rokko::matrix_col_major> locX(dim, 1);
  rokko::localized_matrix<double, rokko::matrix_col_major> locY(dim, 1);
  for (int i = 0; i < dim; ++i) locX(i, 0) = dist(engine);
  for (int i = 0; i < dim; ++i) locY(i, 0) = dist(engine);
  rokko::scatter(locX, vecX, 0);
  rokko::scatter(locY, vecY, 0);

  // local calculation
  double product_local = 0;
  for (int i = 0; i < dim; ++i)
    product_local += locX(i, 0) * locY(i, 0);

  // global calculation
  double product_global = rokko::dot_product(vecX, false, 0, vecY, false, 0);

  int success_local = 1;
  for (int i = 0; i < dim; ++i) {
    if (vecX.is_gindex(i, 0)) {
      std::cerr << i << ' ' << rank << ' ' << product_local << ' ' << product_global << std::endl;
      std::cerr.flush();
      if (std::abs(product_local - product_global) >  10e-12) {
        success_local = 0;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  int success;
  MPI_Allreduce(&success_local, &success, 1, MPI_INT, MPI_PROD, comm);
  ASSERT_TRUE(success);
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
