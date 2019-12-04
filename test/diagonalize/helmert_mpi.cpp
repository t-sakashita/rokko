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
#include <rokko/utility/helmert_matrix.hpp>
#include <rokko/collective.hpp>
#include <rokko/utility/compare_vectors.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

constexpr double eps = 1e-5;

template <typename GRID_MAJOR>
void test(MPI_Comm comm, int dim, std::string const& name, GRID_MAJOR const& grid_major) {
  std::cout << "solver=" << name << std::endl;
  if (name == "eigenexa") {
    std::cout << "warning: skipping test for eigenexa\n";
  } else {
    rokko::parallel_dense_ev solver(name);
    if (solver.is_available_grid_major(grid_major)) {
      solver.initialize(global_argc, global_argv);
      rokko::grid g(comm, grid_major);
      rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(dim, g);
      rokko::distributed_matrix<double, rokko::matrix_col_major> mat(map);
      Eigen::VectorXd diag(dim);
      diag.setLinSpaced(diag.size(), 1, diag.size()); // diag = [1, 2, 3, ..., dim]
      rokko::helmert_matrix::generate_for_given_eigenvalues(mat, diag);
      Eigen::VectorXd w(dim);
      rokko::distributed_matrix<double, rokko::matrix_col_major> Z(map);

      solver.diagonalize(mat, w, Z);

      EXPECT_NEAR(w.sum(), diag.sum(), eps);

      Eigen::MatrixXd locZ(dim,dim);
      constexpr int root = 0;
      rokko::gather(Z, locZ, root);

      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<rokko::matrix_col_major>> u(dim, dim);
      rokko::helmert_matrix::generate(u);
      for (int i = 0; i < dim; ++i) {
        EXPECT_NEAR(diag(i), w(i), eps);

        // The following test utilizes that the fact that all these eigenvectors have freedom of minus sign at most, because all these eigenvalues are simple.
        if (g.get_myrank() == root)
          EXPECT_NEAR(rokko::norm_diff(u.transpose().col(i), locZ.col(i)), 0, eps);
      }

      solver.finalize();
    }
  }
}

TEST(diagonalize, helmert_mpi) {
  MPI_Comm comm = MPI_COMM_WORLD;
  const int dim = 100;

  std::vector<std::string> names;
  if (global_argc == 1) {
    names = rokko::parallel_dense_ev::solvers();
  } else {
    for (int num=1; num < global_argc; ++num) {
      names.push_back(global_argv[num]);
    }
  }

  for(auto name : names) {
    test(comm, dim, name, rokko::grid_col_major);
  }

  for(auto name : names) {
    test(comm, dim, name, rokko::grid_row_major);
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
