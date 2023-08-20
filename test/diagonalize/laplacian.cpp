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
#include <rokko/utility/laplacian_matrix.hpp>
#include <rokko/utility/command_line_parameters.hpp>

#include <gtest/gtest.h>

constexpr double eps = 1e-5;

int global_argc;
char** global_argv;

template<typename MATRIX_MAJOR>
void test(int dim, std::string const& name) {
  rokko::serial_dense_ev solver(name);
  solver.initialize(global_argc, global_argv);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<MATRIX_MAJOR>> mat(dim, dim);
  rokko::laplacian_matrix::generate(mat);
  Eigen::VectorXd eigval(dim);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<MATRIX_MAJOR>> eigvec(dim, dim);

  solver.diagonalize(mat, eigval, eigvec);

  EXPECT_NEAR(eigval.array().inverse().sum(), dim * (dim+1) * 0.5, eps);

  rokko::laplacian_matrix::generate(mat);
  for (int i = 0; i < dim; ++i) {
    const auto w = eigvec.col(i).transpose() * mat * eigvec.col(i);
    EXPECT_NEAR(w, eigval[i], eps);
  }

  solver.finalize();
}

TEST(diagonalize, laplacian) {
  constexpr int dim = 100;
  std::cout << "dimension = " << dim << std::endl;

  const auto names = global_argc == 1 ? rokko::serial_dense_ev::solvers()
    : rokko::get_command_line_args(global_argc, global_argv);

  for(auto const& name : names) {
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
