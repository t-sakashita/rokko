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
#include <rokko/lapack/heev.hpp>
#include <rokko/lapack/lange.hpp>
#include <rokko/eigen3.hpp>
#include <iostream>

int main(int argc, char** argv) {
  constexpr double eps = 1e-5;
  int n = 5;
  if (argc > 1) n = std::stoi(argv[1]);

  // generate matrix
  Eigen::MatrixXcf a = Eigen::MatrixXcf::Random(n, n);
  a += a.adjoint().eval(); // eval() is required to avoid aliasing issue
  std::cout << "Matrix A: " << std::endl << a << std::endl;

  // diagonalization
  Eigen::MatrixXcf u = a;
  Eigen::VectorXf w(n);
  int info = rokko::lapack::heev('V', 'U', u, w);
  if (info) throw std::runtime_error("Error: heev failed");
  std::cout << "Eigenvalues: " << std::endl << w << std::endl;
  std::cout << "Eigenvectors: " << std::endl << u << std::endl;

  // orthogonality check
  Eigen::MatrixXcf check1 = u.adjoint() * u - Eigen::MatrixXcf::Identity(n, n);
  double norm1 = rokko::lapack::lange('F', check1);
  std::cout << "|| U^t U - I || = " << norm1 << std::endl;
  if (norm1 > eps) throw std::runtime_error("Error: orthogonality check");

  // eigenvalue check
  Eigen::MatrixXcf check2 = u.adjoint() * a * u;
  for (int i = 0; i < n; ++i) check2(i, i) -= w(i);
  double norm2 = rokko::lapack::lange('F', check2);
  norm2 = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      norm2 += norm(check2(i, j)) * norm(check2(i, j));
    }
  }
  std::cout << norm2 << std::endl;
  std::cout << "|| U^t A U - diag(w) || = " << norm2 << std::endl;
  if (norm2 > eps) throw std::runtime_error("Error: eigenvalue check");

  return 0;
}
