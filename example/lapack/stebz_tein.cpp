/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/blas.hpp>
#include <rokko/lapack/stebz.hpp>
#include <rokko/lapack/stein.hpp>
#include <rokko/lapack/lange.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/laplacian_matrix.hpp>
#include <iostream>

int main(int argc, char** argv) {
  constexpr double eps = 1e-10;
  int n = 5;
  if (argc > 1) n = std::stoi(argv[1]);

  // generate matrix
  Eigen::VectorXd d(n), e(n-1);
  rokko::laplacian_matrix::generate(d, e);

  // diagonalization
  Eigen::VectorXd w(n);
  constexpr double abstol = 0.;
  int m, nsplit;
  Eigen::VectorXi iblock(n), isplit(n);
  int info = rokko::lapack::stebz('A', 'E', 0., 0., 1, n, abstol, d, e, m, nsplit, w, iblock, isplit);
  Eigen::MatrixXd u = Eigen::MatrixXd::Identity(n,n);
  Eigen::VectorXi ifailv(m);
  info = rokko::lapack::stein(d, e, m, w, u, iblock, isplit, ifailv);
  std::cout << "Eigenvalues: " << std::endl << w << std::endl;
  std::cout << "Eigenvectors: " << std::endl << u << std::endl;

  // orthogonality check
  Eigen::MatrixXd check1 = u.adjoint() * u - Eigen::MatrixXd::Identity(n, n);
  double norm1 = rokko::lapack::lange('F', check1);
  std::cout << "|| U^t U - I || = " << norm1 << std::endl;
  if (norm1 > eps) throw std::runtime_error("Error: orthogonality check");

  // eigenvalue check
  Eigen::MatrixXd a(n,n);
  rokko::laplacian_matrix::generate(a);
  Eigen::MatrixXd check2 = u.adjoint() * a * u;
  check2.diagonal() -= w;
  double norm2 = rokko::lapack::lange('F', check2);
  std::cout << check2 << std::endl;
  std::cout << "|| U^t A U - diag(w) || = " << norm2 << std::endl;
  if (norm2 > eps) throw std::runtime_error("Error: eigenvalue check");

  return 0;
}
