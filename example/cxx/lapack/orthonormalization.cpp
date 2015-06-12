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
  int m = 10;
  int n = 6;
  if (argc > 2) {
    m = boost::lexical_cast<int>(argv[1]);
    n = boost::lexical_cast<int>(argv[2]);
  }
  std::cout << "m = " << m << "\nn = " << n << std::endl;
  int k = std::min(m, n);

  // generate random martix
  matrix_t a = matrix_t::Random(m, n);
  std::cout << "Input random column vectors A:\n" << a << std::endl;

  // orthonormalization
  matrix_t mat = a;
  vector_t tau(k);
  info = LAPACKE_geqrf(LAPACK_COL_MAJOR, m, n, &mat(0, 0), m, &tau(0));
  matrix_t r = mat;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < i; ++j)
      r(i, j) = 0;
  std::cout << "Upper triangle matrix R:\n" << r << std::endl;
  info = LAPACKE_ungqr(LAPACK_COL_MAJOR, m, k, k, &mat(0, 0), m, &tau(0));
  matrix_t q(m, k);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < k; ++j)
      q(i, j) = mat(i, j);
  std::cout << "Orthonormalized column vectors Q:\n" << q << std::endl;

  // check orthonormality
  matrix_t check = q.adjoint() * q;
  std::cout << "Check orthogonality:\n" << check << std::endl;
  std::cout << "Error: | I - Qt Q | = " <<  (matrix_t::Identity(n, n) - check).norm() << std::endl;

  // check decomposition
  std::cout << "Reconstruction of original matrix QR:\n" << q * r << std::endl;
  std::cout << "Error: | A - QR | = " <<  (a - q * r).norm() << std::endl;
}
