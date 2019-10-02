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

#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>

#include <gtest/gtest.h>

TEST(localized_matrix, 123) {

  // M = 1 4 7
  //     2 5 8
  //     3 6 9
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

  int dim = 3;
  rokko::localized_matrix<double> M(dim,dim);
  M << 1,2,3,4,5,6,7,8,9;
  std::cout << M << std::endl;
  double a = 5.0;
  Eigen::VectorXd u(dim);
  u << 1,2,3;
  Eigen::VectorXd v(dim);
  v << 4,5,6;

  Eigen::VectorXd w = a*u+M*v;
  std::cout << w << std::endl;
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
