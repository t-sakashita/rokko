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

#include <rokko/lapack/lange.hpp>
#include <rokko/localized_matrix.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

int main(int argc, char** argv) {
  int m = 5;
  int n = 3;
  
  rokko::slmatrix a = rokko::slmatrix::Random(m, n);
  std::cout << "Matrix A: " << std::endl << a << std::endl;
  std::cout << "|| A || = " << rokko::lapack::lange('F', a) << std::endl;

  rokko::dlmatrix b = rokko::dlmatrix::Random(m, n);
  std::cout << "Matrix B: " << std::endl << b << std::endl;
  std::cout << "|| B || = " << rokko::lapack::lange('F', b) << std::endl;

  rokko::clmatrix c = rokko::clmatrix::Random(m, n);
  std::cout << "Matrix : C" << std::endl << c << std::endl;
  std::cout << "|| C || = " << rokko::lapack::lange('F', c) << std::endl;

  rokko::zlmatrix d = rokko::zlmatrix::Random(m, n);
  std::cout << "Matrix D: " << std::endl << d << std::endl;
  std::cout << "|| D || = " << rokko::lapack::lange('F', d) << std::endl;
  
  return 0;
}
