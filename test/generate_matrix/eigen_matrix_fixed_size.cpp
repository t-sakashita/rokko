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

#include <rokko/eigen3.hpp>

#include <gtest/gtest.h>

TEST(eigen_matrix, fixed_size) {

  // M = 1 2 3
  //     4 5 6
  //     7 8 9
  //
  // a = 5
  //
  // u = 1
  //     2
  //     3
  //
  // v = 4
  //     5
  //     6
  //
  // w = au + Mv = 37
  //               87
  //              137

  constexpr int dim = 3;
  Eigen::Matrix<double,3,3,Eigen::RowMajor> M(dim,dim);
  M << 1,2,3,4,5,6,7,8,9;
  constexpr double a = 5.0;
  Eigen::Vector3d u(dim);
  u << 1,2,3;
  Eigen::Vector3d v(dim);
  v << 4,5,6;

  const Eigen::VectorXd w = a*u+M*v;
  ASSERT_EQ(w[0], 37.);
  ASSERT_EQ(w[1], 87.);
  ASSERT_EQ(w[2], 137.);

  // u^t M = 30 36 42
  ASSERT_EQ( (u.transpose()*M.col(0)), 30.);
  ASSERT_EQ( (u.transpose()*M.col(1)), 36.);
  ASSERT_EQ( (u.transpose()*M.col(2)), 42.);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
