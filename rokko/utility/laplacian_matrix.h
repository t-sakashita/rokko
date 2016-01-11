/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_UTILITY_LAPLACIAN_MATRIX_H
#define ROKKO_UTILITY_LAPLACIAN_MATRIX_H

#include <math.h>

// calculate k-th smallest eigenvalue of dim-dimensional Laplacian matrix (k=0...dim-1)
double rokko_laplacian_matrix_eigenvalue(int dim, int k) {
  return 2 * (1 - cos(M_PI * (2 * k + 1) / (2 * dim + 1)));
}

#endif // ROKKO_UTILITY_LAPLACIAN_MATRIX_H
