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
#include <rokko/lapack.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/localized_matrix.hpp>
#include <boost/lexical_cast.hpp>

int imax(int x, int y) { return (x > y) ? x : y; }

int main(int argc, char** argv) {
  int n = 5;
  if (argc > 1) n = boost::lexical_cast<int>(argv[1]);

  // generate matrix
  rokko::dlmatrix a(n, n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      a(i, j) = n - std::max(i, j);
    }
  }
  std::cout << "Matrix A: " << a.rows() << ' ' << a.cols() << std::endl
            << a << std::endl;

  // diagonalization
  rokko::dlmatrix u = a;
  rokko::dlvector w(n);
  int info = rokko::lapack::syev('V', 'U', u, w);
  std::cout << "Eigenvalues: " << size(w) << std::endl
            << w << std::endl;
  std::cout << "Eigenvectors: " << rows(u) << ' ' << cols(u) << std::endl
            << u << std::endl;

  // orthogonality check
  rokko::dlmatrix t(n, n);
  rokko::blas::gemm(CblasTrans, CblasNoTrans, 1, u, u, 0, t);
  for (int i = 0; i < n; ++i) t(i, i) -= 1;
  double norm2 = 0;
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      norm2 = t(i, j) * t(i, j);
    }
  }
  std::cout << "|| U^t U - I ||^2 = " << norm2 << std::endl;
  if (norm2 > 1e-16) {
    std::cerr << "Error: orthogonality check" << std::endl;
    exit(255);
  }

  // eigenvalue check
  rokko::blas::gemm(CblasNoTrans, CblasNoTrans, 1, a, u, 0, t);
  rokko::blas::gemm(CblasTrans, CblasNoTrans, 1, u, t, 0, a);
  for (int i = 0; i < n; ++i) a(i, i) -= w(i);
  norm2 = 0;
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      norm2 = a(i, j) * a(i, j);
    }
  }
  std::cout << "|| U^t A U - diag(w) ||^2 = " << norm2 << std::endl;
  if (norm2 > 1e-16) {
    fprintf(stderr, "Error: eigenvalue check\n");
    exit(255);
  }

  return 0;
}
