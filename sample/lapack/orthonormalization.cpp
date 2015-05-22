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
#include <rokko/lapacke.hpp>

typedef rokko::localized_vector<double> vector_t;
typedef rokko::localized_matrix<double, rokko::matrix_col_major> matrix_t;

int main(int argc, char *argv[]) {
  int info;
  int m = 10; // m >= n
  int n = 6;
  int k = std::min(m, n);

  // generate random martix
  matrix_t a = matrix_t::Random(m, n);
  std::cout << "Input random column vectors A:\n" << a << std::endl;

  // orthonormalization
  vector_t tau(k);
  info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, m, n, &a(0, 0), m, &tau(0));
  info = LAPACKE_dorgqr(LAPACK_COL_MAJOR, m, n, n, &a(0, 0), m, &tau(0));
  std::cout << "Orthonormalized column vectors V:\n" << a << std::endl;

  // check orthonormality
  matrix_t check = a.adjoint() * a;
  std::cout << "Check orthogonality:\n" << check << std::endl;
  std::cout << "| I - Vt V | = " <<  (matrix_t::Identity(n, n) - check).norm() << std::endl;
}
