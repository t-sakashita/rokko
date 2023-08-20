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

#include <gtest/gtest.h>

template<typename T, typename MATRIX_MAJOR>
void test(int dim) {
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<MATRIX_MAJOR>> mat(dim, dim);
  rokko::laplacian_matrix::generate(mat);
  ASSERT_EQ(mat.trace(), static_cast<T>(2*(dim-1) + 1));  // laplacian matrix consists of integer elements. Hence, its trace is also integer, and no rounding error occurs.
  ASSERT_EQ(mat.sum(), static_cast<T>(1));  // sum of all matrix elements
  ASSERT_TRUE(mat.transpose() == mat);  // checking matrix symmetry
}

TEST(generate_matrix, laplacian_matrix) {
  constexpr int dim = 20;
  std::cout << "dimension = " << dim << std::endl;

  std::cout << "  test for double, row major" << std::endl;
  test<double, rokko::matrix_row_major>(dim);
  std::cout << "  test for double, column major" << std::endl;
  test<double, rokko::matrix_col_major>(dim);
  std::cout << "  test for float, row major" << std::endl;
  test<float, rokko::matrix_row_major>(dim);
  std::cout << "  test for float, column major" << std::endl;
  test<float, rokko::matrix_col_major>(dim);
  std::cout << "  test for int, row major" << std::endl;
  test<int, rokko::matrix_row_major>(dim);
  std::cout << "  test for int, column major" << std::endl;
  test<int, rokko::matrix_col_major>(dim);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
