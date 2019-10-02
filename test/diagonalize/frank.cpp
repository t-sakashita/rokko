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

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

template<typename MATRIX_MAJOR>
void test(int dim, std::string const& name) {
  rokko::serial_dense_ev solver(name);
  solver.initialize(global_argc, global_argv);
  rokko::localized_matrix<double, MATRIX_MAJOR> mat(dim, dim);
  rokko::frank_matrix::generate(mat);
  Eigen::VectorXd eigval(dim);
  rokko::localized_matrix<double, MATRIX_MAJOR> eigvec(dim, dim);

  solver.diagonalize(mat, eigval, eigvec);
  
  double sum = 0;
  for(int i = 0; i < dim; ++i) sum += eigval[i];
  EXPECT_NEAR(sum, dim * (dim+1) * 0.5, 10e-5);
  
  rokko::frank_matrix::generate(mat);
  for (int i = 0; i < dim; ++i) {
    double w = eigvec.col(i).transpose() * mat * eigvec.col(i);
    EXPECT_NEAR(w, eigval[i], 10e-5);
  }

  solver.finalize();
}

TEST(diagonalize, frank) {
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
