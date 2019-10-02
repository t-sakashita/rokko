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

#include <gtest/gtest.h>


template<typename T, typename MATRIX_MAJOR>
void test(int dim) {
  rokko::localized_matrix<T, MATRIX_MAJOR> mat(dim, dim);
  Eigen::VectorXd diag(dim);
  diag.setLinSpaced(diag.size(), 1, diag.size()); // diag = [1, 2, 3, ..., dim]
  rokko::helmert_matrix::generate_for_given_eigenvalues(mat, diag);  
  EXPECT_NEAR(mat.trace(), diag.sum(), 10e-5);
}

TEST(generate_matrix, helmert_matrix) {
  const int dim = 100;
  std::cout << "dimension = " << dim << std::endl;

  std::cout << "  test for row major" << std::endl;
  test<double, rokko::matrix_row_major>(dim);
  std::cout << "  test for column major" << std::endl;
  test<double, rokko::matrix_col_major>(dim);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
