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
#include <rokko/utility/compare_vectors.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

template<typename MATRIX_MAJOR>
void test(int dim, std::string const& name) {
  rokko::serial_dense_ev solver(name);
  solver.initialize(global_argc, global_argv);
  Eigen::VectorXd diag(dim);
  diag.setLinSpaced(diag.size(), 1, diag.size()); // diag = [1, 2, 3, ..., dim]
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<MATRIX_MAJOR>> mat(dim, dim);
  rokko::helmert_matrix::generate_for_given_eigenvalues(mat, diag);
  Eigen::VectorXd eigval(dim);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<MATRIX_MAJOR>> eigvec(dim, dim);

  solver.diagonalize(mat, eigval, eigvec);

  rokko::helmert_matrix::generate_for_given_eigenvalues(mat, diag);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<MATRIX_MAJOR>> u(dim, dim);
  rokko::helmert_matrix::generate(u);

  for (int i = 0; i < dim; ++i) {
    double w = eigvec.col(i).transpose() * mat * eigvec.col(i);
    EXPECT_NEAR(w, eigval(i), 10e-5);

    EXPECT_NEAR(diag(i), eigval(i), 10e-5);

    // The following test utilizes that the fact that all these eigenvectors have freedom of minus sign at most, because all these eigenvalues are simple.
    EXPECT_NEAR(rokko::norm_diff(u.transpose().col(i), eigvec.col(i)), 0, 1e-5);
  }

  solver.finalize();
}

TEST(diagonalize, helmert) {
  const int dim = 100;
  std::cout << "dimension = " << dim << std::endl;

  std::vector<std::string> names;
  if (global_argc == 1) {
    names = rokko::serial_dense_ev::solvers();
  } else {
    for (int num=1; num < global_argc; ++num) {
      names.push_back(global_argv[num]);
    }
  }

  for(auto name : names) {
    std::cout << "solver = " << name << std::endl;
    std::cout << "  test for row major" << std::endl;
    test<rokko::matrix_row_major>(dim, name);
    std::cout << "  test for column major" << std::endl;
    test<rokko::matrix_col_major>(dim, name);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  global_argc = argc;
  global_argv = argv;
  return RUN_ALL_TESTS();
}
