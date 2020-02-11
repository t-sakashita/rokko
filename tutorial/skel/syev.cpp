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
#include <rokko/lapack/syev.hpp>
#include <rokko/lapack/lange.hpp>
#include <rokko/eigen3.hpp>

int main(int argc, char** argv) {
  int n = 5;
  if (argc > 1) n = std::stoi(argv[1]);

  // generate matrix
  Eigen::MatrixXd a(n, n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      a(i, j) = std::min(i, j) + 1;
    }
  }
  std::cout << "Matrix A: " << std::endl << a << std::endl;

  // diagonalization
  Eigen::MatrixXd u = a;
  Eigen::VectorXd w(n);
  int info = rokko::lapack::syev('V', 'U', u, w);
  if (info) throw std::runtime_error("Error: syev failed");
  std::cout << "Eigenvalues: " << std::endl << w << std::endl;
  std::cout << "Eigenvectors: " << std::endl << u << std::endl;

  // orthogonality check
  Eigen::MatrixXd check1 = u.adjoint() * u - Eigen::MatrixXd::Identity(n, n);
  double norm1 = rokko::lapack::lange('F', check1);
  std::cout << "|| U^t U - I || = " << norm1 << std::endl;
  if (norm1 > 1e-10) throw std::runtime_error("Error: orthogonality check");

  // eigenvalue check
  Eigen::MatrixXd check2 = u.adjoint() * a * u;
  for (int i = 0; i < n; ++i) check2(i, i) -= w(i);
  double norm2 = rokko::lapack::lange('F', check2);
  std::cout << "|| U^t A U - diag(w) || = " << norm2 << std::endl;
  if (norm2 > 1e-10) throw std::runtime_error("Error: eigenvalue check");

  return 0;
}
