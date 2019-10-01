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

#include <rokko/grid.hpp>
#include <rokko/solver.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/utility/frank_matrix.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

TEST(generate_matrix, frank_mpi) {
  unsigned int dim = 1000;
  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid g(comm);
  for(auto name : rokko::parallel_dense_ev::solvers()) {
    rokko::parallel_dense_ev solver(name);
    solver.initialize(global_argc, global_argv);
    rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(dim, g);
    rokko::distributed_matrix<double,rokko::matrix_col_major> mat(map);
    rokko::frank_matrix::generate(mat);
    solver.finalize();
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
