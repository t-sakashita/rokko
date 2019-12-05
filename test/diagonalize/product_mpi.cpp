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

#include <rokko/rokko.hpp>
#include <rokko/utility/frank_matrix.hpp>
#include <boost/lexical_cast.hpp>

#include <gtest/gtest.h>

constexpr double eps = 1e-11;

int global_argc;
char** global_argv;

TEST(product_mpi, product_mpi) {
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  int dim = 100;
  if (global_argc > 1) {
    dim = boost::lexical_cast<int>(global_argv[1]);
  }

  if (rank == 0) std::cout << "dimension = " << dim << std::endl;
  rokko::parallel_dense_ev solver(rokko::parallel_dense_ev::default_solver());
  rokko::grid g(comm);
  rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(dim, g);
  rokko::distributed_matrix<double, rokko::matrix_col_major> matA(map);
  rokko::distributed_matrix<double, rokko::matrix_col_major> matB(map);
  rokko::distributed_matrix<double, rokko::matrix_col_major> matC(map);
  rokko::frank_matrix::generate(matA);
  rokko::frank_matrix::generate(matB);
  rokko::product(1.0, matA, false, matB, false, 0, matC);
  matC.print();
  // calculate trace
  double sum_local = 0;
  for (int i = 0; i < dim; ++i) {
    if (matC.is_gindex(i, i)) sum_local += matC.get_global(i, i);
  }
  double sum_global = 0;
  MPI_Allreduce(&sum_local, &sum_global, 1, MPI_DOUBLE, MPI_SUM, comm);
  if (rank == 0) std::cout << "trace of distributed matrix = " << sum_global << std::endl;

  Eigen::MatrixXd lmatA(dim, dim);
  rokko::frank_matrix::generate(lmatA);
  Eigen::MatrixXd lmatC = lmatA * lmatA;
  if (rank == 0) std::cout << lmatC << std::endl;
  double sum = lmatC.trace();
  if (rank == 0) std::cout << "trace of eigen matrix = " << sum << std::endl;

  if (rank == 0) EXPECT_NEAR(sum_global, sum, eps);
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
