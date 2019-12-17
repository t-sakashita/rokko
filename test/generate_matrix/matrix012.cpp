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
#include <rokko/utility/matrix012.hpp>

#include <gtest/gtest.h>

template <typename T, int MAJOR=Eigen::ColMajor>
using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1, MAJOR>;

// for general rectangular matrix
template<typename T, typename MATRIX_MAJOR>
void test(int dim1, int dim2) {
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<MATRIX_MAJOR>> mat(dim1, dim2);
  rokko::matrix012::generate(mat);
  Eigen::Map<VectorX<T>> vec(mat.data(), mat.size());  // flattening mat

  ASSERT_TRUE(vec == VectorX<T>::LinSpaced(vec.size(), 0, vec.size()-1));
}

// for square matrix
template<typename T, typename MATRIX_MAJOR>
void test(int dim) {
  test<T,MATRIX_MAJOR>(dim, dim);
}

template<typename T>
void test_type() {
  test<T, rokko::matrix_row_major>(10);
  test<T, rokko::matrix_col_major>(10);
  test<T, rokko::matrix_row_major>(3, 10);
  test<T, rokko::matrix_col_major>(3, 10);
  test<T, rokko::matrix_row_major>(10, 3);
  test<T, rokko::matrix_col_major>(10, 3);
}

TEST(generate_matrix, minij_matrix) {
  // for int
  std::cout << "int" << std::endl;
  test_type<int>();
  // for double
  std::cout << "dobule" << std::endl;
  test_type<double>();
  // for flaot
  std::cout << "float" << std::endl;
  test_type<float>();
}

int main(int argc, char** argv) {
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  result = RUN_ALL_TESTS();
  return result;
}
