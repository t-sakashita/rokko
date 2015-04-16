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

#include <boost/random.hpp>
#include <rokko/rokko.hpp>
#include <lapacke.h>

int main(int argc, char *argv[]) {
  int info;
  const int dim = 6;

  // generate random martix
  rokko::localized_matrix<> mat(dim, dim);
  boost::mt19937 eng(1234ul);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
    rng(eng, boost::normal_distribution<>());
  for (int j = 0; j < dim; ++j)
    for (int i = 0; i < dim; ++i)
      mat(i, j) = rng();
  std::cout << "Input random column vectors:\n" << mat << std::endl;

  // orthogonalization
  rokko::localized_vector tau(dim);
  info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, dim, dim, &mat(0, 0), dim, &tau(0));
  info = LAPACKE_dorgqr(LAPACK_COL_MAJOR, dim, dim, dim, &mat(0, 0), dim, &tau(0));
  std::cout << "Orthonormalized column vectors:\n" << mat << std::endl;

  // check orthogonality
  std::cout << "Check orthogonality:\n" << mat.transpose() * mat << std::endl;
}
