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

TEST(GemvTest, GemvTest) {
  std::size_t n = 16;
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(n, n);
  Eigen::VectorXd x = Eigen::VectorXd::Random(n);
  Eigen::VectorXd y = Eigen::VectorXd::Random(n);
  Eigen::VectorXd yr = y;
  double alpha = 2.3;
  double beta = 0.5;
  rokko::blas::gemv(CblasNoTrans, alpha, a, x, 1, beta, yr, 1);

  Eigen::VectorXd yc(n);
  for (std::size_t i = 0; i < n; ++i) {
    yc(i) = beta * y(i);
    for (std::size_t j = 0; j < n; ++j) {
      yc(i) += alpha * a(i, j) * x(j);
    }
    EXPECT_NEAR(yc(i), yr(i), 1.0e-12);
  }
}
