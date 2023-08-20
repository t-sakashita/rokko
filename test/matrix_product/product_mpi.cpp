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
#include <rokko/collective.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

TEST(product_mpi, product_mpi) {
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  const int dim = (global_argc > 1) ? std::stoi(global_argv[1]) : 100;

  if (rank == 0) std::cout << "dimension = " << dim << std::endl;
  rokko::parallel_dense_ev solver(rokko::parallel_dense_ev::default_solver());
  rokko::grid g(comm);
  rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(dim, g);
  rokko::distributed_matrix<double, rokko::matrix_col_major> matA(map);
  rokko::distributed_matrix<double, rokko::matrix_col_major> matC(map);
  rokko::frank_matrix::generate(matA);
  rokko::product(1.0, matA, false, matA, false, 0, matC);
  matC.print();
  double trace_distributed = rokko::trace(matC);
  const int dim_proc = (rank == 0) ? dim : 0;
  Eigen::MatrixXd g_matC(dim_proc, dim_proc);
  rokko::gather(matC, g_matC, 0);

  if (rank == 0) {
    std::cout << "trace of distributed matrix = " << trace_distributed << std::endl;
    Eigen::MatrixXd lmatA(dim, dim);
    rokko::frank_matrix::generate(lmatA);
    const Eigen::MatrixXd lmatC = lmatA * lmatA;
    std::cout << lmatC << std::endl;
    const double trace = lmatC.trace();
    std::cout << "trace of eigen matrix = " << trace << std::endl;
    ASSERT_EQ(trace_distributed, trace);
    ASSERT_TRUE(g_matC == lmatC);
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
