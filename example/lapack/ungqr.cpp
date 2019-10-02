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
#include <rokko/lapack/geqrf.hpp>
#include <rokko/lapack/ungqr.hpp>
#include <rokko/lapack/lange.hpp>
#include <rokko/localized_matrix.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
  int m = 10;
  int n = 6;
  if (argc > 2) {
    m = boost::lexical_cast<int>(argv[1]);
    n = boost::lexical_cast<int>(argv[2]);
  }
  int k = std::min(m, n);
  if (m < n) throw std::invalid_argument("Error: m < n");
  
  // generate matrix
  Eigen::VectorXcd a = rokko::zlmatrix::Random(m, n);
  std::cout << "Matrix A: " << std::endl << a << std::endl;

  // orthonormaliation
  rokko::zlmatrix mat = a;
  Eigen::VectorXcd tau(k);
  int info = rokko::lapack::geqrf(mat, tau);
  if (info) throw std::runtime_error("Error: geqrf failed");

  Eigen::VectorXcd r = mat;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < i; ++j)
      r(i, j) = 0;
  std::cout << "Upper triangle matrix R:" << std::endl << r << std::endl;
  info = rokko::lapack::orgqr(k, mat, tau);
  rokko::zlmatrix q(m, k);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < k; ++j)
      q(i, j) = mat(i, j);
  std::cout << "Orthonormalized column vectors Q:" << std::endl << q << std::endl;
  
  // orthogonality check
  rokko::zlmatrix check1 = q.adjoint() * q - rokko::zlmatrix::Identity(k, k);
  double norm1 = rokko::lapack::lange('F', check1);
  std::cout << "|| Q^t Q - I || = " << norm1 << std::endl;
  if (norm1 > 1e-10) throw std::runtime_error("Error: orthogonality check");

  // solution check
  rokko::zlmatrix check2 = q * r - a;
  double norm2 = rokko::lapack::lange('F', check2);
  std::cout << "|| Q R - A || = " << norm2 << std::endl;
  if (norm2 > 1e-10) throw std::runtime_error("Error: solution check");

  return 0;
}
