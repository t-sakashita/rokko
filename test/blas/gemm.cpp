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

TEST(GemmTest, GemmTest) {
  std::size_t n = 16;
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(n, n);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(n, n);
  Eigen::MatrixXd c = Eigen::MatrixXd::Random(n, n);
  Eigen::MatrixXd cr = c;
  double alpha = 2.1;
  double beta = 0.8;

  rokko::blas::gemm(CblasNoTrans, CblasNoTrans, alpha, a, b, beta, cr);

  Eigen::MatrixXd cc(n, n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      cc(i, j) = beta * c(i, j);
      for (std::size_t k = 0; k < n; ++k) {
        cc(i, j) += alpha * a(i, k) * b(k, j);
      }
      EXPECT_NEAR(cc(i, j), cr(i, j), eps);
    }
  }
}
