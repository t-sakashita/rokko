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
#include <rokko/utility/minij_matrix.hpp>

#include <gtest/gtest.h>

constexpr double eps = 1e-5;

int global_argc;
char** global_argv;

TEST(diagonalize, minij_mpi) {
  MPI_Comm comm = MPI_COMM_WORLD;
  constexpr int dim = 10;

  std::vector<std::string> names;
  if (global_argc == 1) {
    names = rokko::parallel_dense_ev::solvers();
  } else {
    for (int num=1; num < global_argc; ++num) {
      names.emplace_back(global_argv[num]);
    }
  }

  for(auto name : names) {
    std::cout << "solver=" << name << std::endl;
    rokko::parallel_dense_ev solver(name);
    if (solver.is_available_grid_major(rokko::grid_col_major)) {
      solver.initialize(global_argc, global_argv);
      rokko::grid g(comm, rokko::grid_col_major);
      rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(dim, g);
      rokko::distributed_matrix<double, rokko::matrix_col_major> mat(map);
      rokko::minij_matrix::generate(mat);
      Eigen::VectorXd w(dim);
      rokko::distributed_matrix<double, rokko::matrix_col_major> Z(map);

      solver.diagonalize(mat, w, Z);

      EXPECT_NEAR(w.sum(), dim * (dim+1) * 0.5, eps);
      
      solver.finalize();
    }
  }


  for(auto name : names) {
    std::cout << "solver=" << name << std::endl;
    rokko::parallel_dense_ev solver(name);
    if (solver.is_available_grid_major(rokko::grid_row_major)) {
      solver.initialize(global_argc, global_argv);
      rokko::grid g(comm, rokko::grid_row_major);
      rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(dim, g);
      rokko::distributed_matrix<double, rokko::matrix_col_major> mat(map);
      rokko::minij_matrix::generate(mat);
      Eigen::VectorXd w(dim);
      rokko::distributed_matrix<double, rokko::matrix_col_major> Z(map);
      
      solver.diagonalize(mat, w, Z);

      EXPECT_NEAR(w.sum(), dim * (dim+1) * 0.5, eps);

      solver.finalize();
    }
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
