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
#include <rokko/lapack.hpp>
#include <boost/lexical_cast.hpp>

typedef rokko::localized_vector<double> vector_t;
typedef rokko::localized_matrix<double, rokko::matrix_col_major> matrix_t;

int main(int argc, char *argv[]) {
  int info;
  int m = 6;
  int n = 4;
  if (argc > 2) {
    m = boost::lexical_cast<int>(argv[1]);
    n = boost::lexical_cast<int>(argv[2]);
  }
  std::cout << "m = " << m << "\nn = " << n << std::endl;
  int k = std::min(m, n);

  // generate random martix
  matrix_t mat = matrix_t::Random(m, n);
  std::cout << "Input random matrix A:\n" << mat << std::endl;

  // singular value decomposition
  matrix_t a = mat; // 'a' will be destroyed by dgesvd
  vector_t s(k), superb(k-1);
  matrix_t u(m, k), vt(k, n);
  info = rokko::lapack::gesvd('S', 'S', a, s, u, vt, superb);
  std::cout << "U:\n" << u << std::endl;
  std::cout << "S:\n" << s << std::endl;
  std::cout << "Vt:\n" << vt << std::endl;

  // check correctness of SVD
  matrix_t smat = matrix_t::Zero(k, k);
  for (int i = 0; i < k; ++i) smat(i, i) = s(i);
  matrix_t check = u * smat * vt;
  std::cout << "U * S * Vt:\n" << check << std::endl;
  std::cout << "| A - U * S * Vt | = " << (mat - check).norm() << std::endl;
}
