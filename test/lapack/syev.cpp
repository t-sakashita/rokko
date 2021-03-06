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
constexpr double eps_float = 1e-4;

TEST(syev, dsyev) {
  int n = 32;

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
  EXPECT_NEAR(0.0, norm1, eps);

  // eigenvalue check
  Eigen::MatrixXd check2 = u.adjoint() * a * u;
  for (int i = 0; i < n; ++i) check2(i, i) -= w(i);
  double norm2 = rokko::lapack::lange('F', check2);
  EXPECT_NEAR(0.0, norm2, eps);
}

TEST(syev, ssyev) {
  int n = 32;

  // generate matrix
  Eigen::MatrixXf a = Eigen::MatrixXf::Random(n, n);
  a += a.transpose().eval();
  
  // diagonalization
  Eigen::MatrixXf u = a;
  Eigen::VectorXf w(n);
  int info = rokko::lapack::syev('V', 'U', u, w);
  EXPECT_EQ(0, info);

  // orthogonality check
  Eigen::MatrixXf check1 = u.adjoint() * u - Eigen::MatrixXf::Identity(n, n);
  float norm1 = rokko::lapack::lange('F', check1);
  EXPECT_NEAR(0.0, norm1, eps_float);

  // eigenvalue check
  Eigen::MatrixXf check2 = u.adjoint() * a * u;
  for (int i = 0; i < n; ++i) check2(i, i) -= w(i);
  float norm2 = rokko::lapack::lange('F', check2);
  EXPECT_NEAR(0.0, norm2, eps_float);
}
