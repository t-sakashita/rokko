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

#include <rokko/lapack.hpp>
#include <rokko/localized_matrix.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

int main(int argc, char** argv) {
  int m = 3;
  int n = 5;
  
  // generate matrix
  rokko::zlmatrix a(m, n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < m; ++i) {
      a(i, j) = std::complex<double>(i, j);
    }
  }
  std::cout << "Matrix A: " << rows(a) << ' ' << cols(a) << std::endl
            << a << std::endl;

  double norm = rokko::lapack::lange('F', a);
  std::cout << "|| A || = " << norm << std::endl;

  return 0;
}
