/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <gtest/gtest.h>
#include <rokko/eigen3.hpp>
#include <rokko/blas.hpp>

constexpr double eps = 1e-12;

TEST(GemvTest, GemvTest) {
  constexpr std::size_t n = 16;
  const Eigen::MatrixXd a = Eigen::MatrixXd::Random(n, n);
  const Eigen::VectorXd x = Eigen::VectorXd::Random(n);
  const Eigen::VectorXd y = Eigen::VectorXd::Random(n);
  Eigen::VectorXd yr = y;
  constexpr double alpha = 2.3;
  constexpr double beta = 0.5;
  rokko::blas::gemv(CblasNoTrans, alpha, a, x, 1, beta, yr, 1);

  Eigen::VectorXd yc(n);
  for (std::size_t i = 0; i < n; ++i) {
    yc(i) = beta * y(i);
    for (std::size_t j = 0; j < n; ++j) {
      yc(i) += alpha * a(i, j) * x(j);
    }
    EXPECT_NEAR(yc(i), yr(i), eps);
  }
}
