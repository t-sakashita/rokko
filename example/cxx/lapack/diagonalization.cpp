/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <rokko/lapack.h>
#include <boost/lexical_cast.hpp>

typedef rokko::localized_vector<double> vector_t;
typedef rokko::localized_matrix<double, rokko::matrix_col_major> matrix_t;

int main(int argc, char *argv[]) {
  int info;
  int n = 6;
  if (argc > 1) n = boost::lexical_cast<int>(argv[1]);
  std::cout << "n = " << n << std::endl;

  // generate symmmetric random martix
  matrix_t mat = matrix_t::Random(n, n);
  mat += mat.adjoint().eval(); // eval() is required to avoid aliasing issue
  std::cout << "Input random matrix A:\n" << mat << std::endl;

  // diagonalization
  vector_t w(n);
  matrix_t v = mat; // v will be overwritten by eigenvectors
  info = LAPACKE_heev(LAPACK_COL_MAJOR, 'V', 'U', n, &v(0, 0), n, &w(0));
  std::cout << "Eigenvalues:\n" << w << std::endl;
  std::cout << "Eigenvectors:\n" << v << std::endl;

  // check correctness of diagonalization
  matrix_t wmat = matrix_t::Zero(n, n);
  for (int i = 0; i < n; ++i) wmat(i, i) = w(i);
  matrix_t check = v.adjoint() * mat * v;
  std::cout << "Vt * A * V:\n" << check << std::endl;
  std::cout << "| W - Vt * A * V | = " << (wmat - check).norm() << std::endl;
}
