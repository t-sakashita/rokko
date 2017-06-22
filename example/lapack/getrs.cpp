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

  // generate matrix and rhs vector
  rokko::dlmatrix a(n, n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      a(i, j) = n - std::max(i, j);
    }
  }
  std::cout << "Matrix A: " << rows(a) << ' ' << cols(a) << std::endl
            << a << std::endl;
  rokko::dlvector b(n);
  for (int i = 0; i < n; ++i) b(i) = i * i + 1;
  std::cout << "Vector b: " << size(b) << std::endl
            << b << std::endl;

  // solve linear equation
  rokko::dlmatrix lu = a;
  rokko::dlvector x = b;
  rokko::ilvector ipiv(n);
  int info = rokko::lapack::getrf(lu, ipiv);
  info = rokko::lapack::getrs('n', 1, lu, ipiv, x);
  std::cout << "Solution x: " << size(x) << std::endl
            << x << std::endl;

  /* solution check */
  rokko::blas::gemv(CblasNoTrans, 1, a, x, 1, -1, b, 1);
  double norm2 = 0;
  for (int i = 0; i < n; ++i) norm2 = b(i) * b(i);
  std::cout << "|| A x - b ||^2 = " << norm2 << std::endl;
  if (norm2 > 1e-16) {
    std::cerr << "Error: solution check" << std::endl;
    exit(255);
  }

  return 0;
}
