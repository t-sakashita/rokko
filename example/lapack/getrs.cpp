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

#include <rokko/blas.hpp>
#include <rokko/lapack/getrf.hpp>
#include <rokko/lapack/getrs.hpp>
#include <rokko/eigen3.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

int main(int argc, char** argv) {
  constexpr double eps = 1e-10;
  int n = 5;
  if (argc > 1) n = boost::lexical_cast<int>(argv[1]);

  // generate matrix and rhs vector
  Eigen::MatrixXd a(n, n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      a(i, j) = std::min(i, j) + 1;
    }
  }
  std::cout << "Matrix A: " << std::endl << a << std::endl;
  Eigen::VectorXd b(n);
  for (int i = 0; i < n; ++i) b(i) = i * i + 1;
  std::cout << "Vector b: " << std::endl << b << std::endl;

  // solve linear equation
  Eigen::MatrixXd lu = a;
  Eigen::VectorXd x = b;
  Eigen::VectorXi ipiv(n);
  rokko::lapack::getrf(lu, ipiv);
  rokko::lapack::getrs('n', 1, lu, ipiv, x);
  std::cout << "Solution x: " << std::endl << x << std::endl;

  /* solution check */
  Eigen::VectorXd check = a * x - b;
  double norm = check.norm();
  std::cout << "|| A x - b || = " << norm << std::endl;
  if (norm > eps) throw std::runtime_error("Error: solution check");

  return 0;
}
