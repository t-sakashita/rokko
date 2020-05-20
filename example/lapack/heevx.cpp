/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/blas.hpp>
#include <rokko/lapack/syevx.hpp>
#include <rokko/lapack/lange.hpp>
#include <rokko/eigen3.hpp>
#include <iostream>

int main(int argc, char** argv) {
  constexpr double eps = 1e-10;
  int n = 5;
  if (argc > 1) n = std::stoi(argv[1]);

  // generate matrix
  Eigen::MatrixXcd a(n, n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      a(i, j) = std::min(i, j) + 1;
    }
  }
  Eigen::MatrixXcd a0 = a;
  std::cout << "Matrix A: " << std::endl << a << std::endl;

  // diagonalization
  Eigen::VectorXd w(n);
  Eigen::MatrixXcd u(n,n);
  constexpr double vl = 0, vu = 0;
  constexpr int il = 0, iu = 0;
  constexpr double abstol = 0.;
  int m;
  Eigen::VectorXi ifail(n);

  int info = rokko::lapack::heevx('V', 'A', 'U', a, vl, vu, il, iu, abstol, m, w, u, ifail);
  if (info) throw std::runtime_error("Error: heevx failed");
  std::cout << "Eigenvalues: " << std::endl << w << std::endl;
  std::cout << "Eigenvectors: " << std::endl << u << std::endl;

  // orthogonality check
  Eigen::MatrixXcd check1 = u.adjoint() * u - Eigen::MatrixXcd::Identity(n, n);
  double norm1 = rokko::lapack::lange('F', check1);
  std::cout << "|| U^t U - I || = " << norm1 << std::endl;
  if (norm1 > eps) throw std::runtime_error("Error: orthogonality check");

  // eigenvalue check
  Eigen::MatrixXcd check2 = u.adjoint() * a0 * u;
  check2.diagonal() -= w;
  double norm2 = rokko::lapack::lange('F', check2);
  std::cout << "|| U^t A U - diag(w) || = " << norm2 << std::endl;
  if (norm2 > eps) throw std::runtime_error("Error: eigenvalue check");

  return 0;
}
