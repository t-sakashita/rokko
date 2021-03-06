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

constexpr double eps = 1e-10;

TEST(SyevTest, SyevTest) {
  int n = 10;

  // generate matrix
  Eigen::MatrixXcd a = Eigen::MatrixXcd::Random(n, n);
  a += a.adjoint().eval();

  // diagonalization
  Eigen::MatrixXcd u = a;
  Eigen::VectorXd w(n);
  int info = rokko::lapack::heev('V', 'U', u, w);
  EXPECT_EQ(0, info);

  // orthogonality check
  Eigen::MatrixXcd check1 = u.adjoint() * u - Eigen::MatrixXcd::Identity(n, n);
  double norm1 = rokko::lapack::lange('F', check1);
  EXPECT_NEAR(0.0, norm1, eps);

  // eigenvalue check
  Eigen::MatrixXcd check2 = u.adjoint() * a * u;
  for (int i = 0; i < n; ++i) check2(i, i) -= w(i);
  double norm2 = rokko::lapack::lange('F', check2);
  EXPECT_NEAR(0.0, norm2, eps);
}
