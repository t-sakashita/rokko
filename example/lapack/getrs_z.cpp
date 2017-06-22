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

int main(int argc, char** argv) {
  int n = 5;
  if (argc > 1) n = boost::lexical_cast<int>(argv[1]);

  // generate matrix and rhs vector
  rokko::zlmatrix a = rokko::zlmatrix::Random(n, n);
  std::cout << "Matrix A: " << std::endl << a << std::endl;
  rokko::zlvector b(n);
  for (int i = 0; i < n; ++i) b(i) = i * i + 1;
  std::cout << "Vector b: " << std::endl << b << std::endl;

  // solve linear equation
  rokko::zlmatrix lu = a;
  rokko::zlvector x = b;
  rokko::ilvector ipiv(n);
  int info = rokko::lapack::getrf(lu, ipiv);
  info = rokko::lapack::getrs('n', 1, lu, ipiv, x);
  std::cout << "Solution x: " << std::endl << x << std::endl;

  /* solution check */
  rokko::zlvector check = a * x - b;
  double norm = check.norm();
  std::cout << "|| A x - b || = " << norm << std::endl;
  if (norm > 1e-10) throw std::runtime_error("Error: solution check");

  return 0;
}