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
#include <rokko/lapack/sterf.hpp>
#include <rokko/lapack/lange.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/laplacian_matrix.hpp>

int main(int argc, char** argv) {
  constexpr double eps = 1e-10;
  int n = 5;
  if (argc > 1) n = std::stoi(argv[1]);

  // generate matrix
  Eigen::VectorXd d(n), e(n-1);
  rokko::laplacian_matrix::generate(d, e);

  // diagonalization
  Eigen::VectorXd w = d;
  int info = rokko::lapack::sterf(w, e);
  std::cout << "Eigenvalues: " << std::endl << w << std::endl;

  return 0;
}
