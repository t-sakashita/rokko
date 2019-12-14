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
#include <rokko/utility/minij_matrix.hpp>
#include <rokko/collective.hpp>

#include <gtest/gtest.h>

constexpr double eps = 1e-5;

int global_argc;
char** global_argv;

TEST(generate_matrix, minij_mpi) {
  constexpr int dim = 1000;
  constexpr int root_proc = 0;
  rokko::grid g(MPI_COMM_WORLD);

  for(auto const& name : rokko::parallel_dense_ev::solvers()) {
    rokko::parallel_dense_ev solver(name);
    solver.initialize(global_argc, global_argv);
    rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(dim, g);
    rokko::distributed_matrix<double,rokko::matrix_col_major> mat(map);
    rokko::minij_matrix::generate(mat);
    double trace_mat = trace(mat);

    if (g.get_myrank() == root_proc) {
      constexpr double trace_th = dim * (dim+1) * 0.5;
      std::cout << "trace = " << trace_mat << " sum_eigenvalues=" << trace_th << std::endl;
      EXPECT_NEAR(trace_mat, trace_th, trace_th * eps);  // expect relative error < eps
    }

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
