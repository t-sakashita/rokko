/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/eigen3.hpp>
#include <type_traits>

#include <gtest/gtest.h>

TEST(matrix, major_template_transformation) {
  static_assert( std::is_same_v<rokko::matrix_major<Eigen::MatrixXd>, rokko::matrix_col_major>, "rokko::matrix_major<Eigen::MatrixXd> != rokko::matrix_col_major");
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
