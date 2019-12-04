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
#include <cmath>
#include <rokko/eigen3.hpp>
#include <rokko/lapack.hpp>

constexpr double eps = 1e-10;
constexpr double eps_float = 1e-5;

TEST(lange, dlange) {
  int m = 5;
  int n = 3;

  // generate matrix
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(m, n);
  
  std::cout << a << std::endl;
  
  double norm, expect;
  
  norm = rokko::lapack::lange('M', a);
  expect = 0.0;
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i)
      expect = std::max(expect, std::abs(a(i, j)));
  std::cout << "max norm of A = " << norm << " (expect: " << expect << ")" << std::endl;
  EXPECT_NEAR(expect, norm, eps);

  norm = rokko::lapack::lange('1', a);
  expect = 0.0;
  for (int j = 0; j < n; ++j) {
    double sum = 0.0;
    for (int i = 0; i < m; ++i) sum += std::abs(a(i, j));
    expect = std::max(expect, sum);
  }
  std::cout << "1-norm of A = " << norm << " (expect: " << expect << ")" << std::endl;
  EXPECT_NEAR(expect, norm, eps);

  norm = rokko::lapack::lange('I', a);
  expect = 0.0;
  for (int i = 0; i < m; ++i) {
    double sum = 0.0;
    for (int j = 0; j < n; ++j) sum += std::abs(a(i, j));
    expect = std::max(expect, std::abs(sum));
  }
  std::cout << "infinity norm of A = " << norm << " (expect: " << expect << ")" << std::endl;
  EXPECT_NEAR(expect, norm, eps);

  norm = rokko::lapack::lange('F', a);
  expect = 0.0;
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i) expect += a(i, j) * a(i, j);
  expect = std::sqrt(expect);
  std::cout << "Frobenius norm of A = " << norm << " (expect: " << expect << ")" << std::endl;
  EXPECT_NEAR(expect, norm, eps);
}

TEST(lange, slange) {
  int m = 5;
  int n = 3;

  // generate matrix
  Eigen::MatrixXf a = Eigen::MatrixXf::Random(m, n);
  
  std::cout << a << std::endl;
  
  float norm, expect;
  
  norm = rokko::lapack::lange('M', a);
  expect = 0.0;
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i)
      expect = std::max(expect, std::abs(a(i, j)));
  std::cout << "max norm of A = " << norm << " (expect: " << expect << ")" << std::endl;
  EXPECT_NEAR(expect, norm, eps_float);

  norm = rokko::lapack::lange('1', a);
  expect = 0.0;
  for (int j = 0; j < n; ++j) {
    float sum = 0.0;
    for (int i = 0; i < m; ++i) sum += std::abs(a(i, j));
    expect = std::max(expect, sum);
  }
  std::cout << "1-norm of A = " << norm << " (expect: " << expect << ")" << std::endl;
  EXPECT_NEAR(expect, norm, eps_float);

  norm = rokko::lapack::lange('I', a);
  expect = 0.0;
  for (int i = 0; i < m; ++i) {
    float sum = 0.0;
    for (int j = 0; j < n; ++j) sum += std::abs(a(i, j));
    expect = std::max(expect, std::abs(sum));
  }
  std::cout << "infinity norm of A = " << norm << " (expect: " << expect << ")" << std::endl;
  EXPECT_NEAR(expect, norm, eps_float);

  norm = rokko::lapack::lange('F', a);
  expect = 0.0;
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i) expect += a(i, j) * a(i, j);
  expect = std::sqrt(expect);
  std::cout << "Frobenius norm of A = " << norm << " (expect: " << expect << ")" << std::endl;
  EXPECT_NEAR(expect, norm, eps_float);
}
