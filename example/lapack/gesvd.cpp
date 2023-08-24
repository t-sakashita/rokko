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
#include <rokko/lapack/gesvd.hpp>
#include <rokko/lapack/lange.hpp>
#include <rokko/eigen3.hpp>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
  constexpr double eps = 1e-10;
  int m = 3;
  int n = 5;
  if (argc > 2) {
    m = std::stoi(argv[1]);
    n = std::stoi(argv[2]);
  }
  int r = std::min(m, n);

  // generate matrix
  Eigen::MatrixXd a(m, n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < m; ++i) {
      a(i, j) = n - std::max(i, j);
    }
  }
  std::cout << "Matrix A: " << std::endl << a << std::endl;

  // singular value decomposition
  Eigen::VectorXd s(r);
  Eigen::MatrixXd u(m, r), vt(r, n);
  Eigen::MatrixXd t = a;
  Eigen::VectorXd superb(r-1);
  int info = rokko::lapack::gesvd('S', 'S', t, s, u, vt, superb);
  if (info) throw std::runtime_error("Error: gesvd failed");
  std::cout << "Matrix U: " << std::endl << u << std::endl;
  std::cout << "Singular values: " << std::endl << s << std::endl;
  std::cout << "Matrix Vt: " << std::endl << vt << std::endl;

  // orthogonality check
  Eigen::MatrixXd check1 = u.adjoint() * u - Eigen::MatrixXd::Identity(r, r);
  double norm1 = rokko::lapack::lange('F', check1);
  std::cout << "|| U^t U - I || = " << norm1 << std::endl;
  if (norm1 > eps) throw std::runtime_error("Error: orthogonality check");

  Eigen::MatrixXd check2 = vt * vt.adjoint() - Eigen::MatrixXd::Identity(r, r);
  double norm2 = rokko::lapack::lange('F', check2);
  std::cout << "|| V^t V - I || = " << norm2 << std::endl;
  if (norm2 > eps) throw std::runtime_error("Error: orthogonality check");

  // solution check
  Eigen::MatrixXd check3 = u.adjoint() * a * vt.adjoint();
  for (int i = 0; i < r; ++i) check3(i, i) -= s(i);
  double norm3 = rokko::lapack::lange('F', check3);
  std::cout << "|| U^t A V - diag(S) || = " << norm3 << std::endl;
  if (norm3 > eps) throw std::runtime_error("Error: solution check");

  return 0;
}
