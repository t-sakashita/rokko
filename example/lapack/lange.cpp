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
  
  rokko::dlmatrix a = rokko::dlmatrix::Random(m, n);
  std::cout << "Matrix A: " << std::endl << a << std::endl;
  std::cout << "|| A || = " << rokko::lapack::lange('F', a) << std::endl;

  rokko::zlmatrix b = rokko::zlmatrix::Random(m, n);
  std::cout << "Matrix B: " << std::endl << b << std::endl;
  std::cout << "|| B || = " << rokko::lapack::lange('F', b) << std::endl;

  return 0;
}
