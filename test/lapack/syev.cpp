/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2017 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <gtest/gtest.h>
#include <rokko/eigen3.hpp>
#include <rokko/lapack.hpp>

TEST(SyevTest, SyevTest) {
  int n = 10;

  // generate matrix
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(n, n);
  a += a.transpose().eval();
  
  // diagonalization
  Eigen::MatrixXd u = a;
  Eigen::VectorXd w(n);
  int info = rokko::lapack::syev('V', 'U', u, w);
  EXPECT_EQ(0, info);

  // orthogonality check
  Eigen::MatrixXd check1 = u.adjoint() * u - Eigen::MatrixXd::Identity(n, n);
  double norm1 = rokko::lapack::lange('F', check1);
  EXPECT_NEAR(0.0, norm1, 1e-10);

  // eigenvalue check
  Eigen::MatrixXd check2 = u.adjoint() * a * u;
  for (int i = 0; i < n; ++i) check2(i, i) -= w(i);
  double norm2 = rokko::lapack::lange('F', check2);
  EXPECT_NEAR(0.0, norm2, 1e-10);
}
