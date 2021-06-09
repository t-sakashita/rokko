/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <rokko/cpblas.h>
#include <rokko/lapack/complex_cast.hpp>

namespace rokko {
namespace pblas {

using rokko::lapack::complex_cast;

inline void pcopy(int n, const float* X, int ix, int jx, const int* descX, int incX,
                  float* Y, int iy, int jy, const int* descY, int incY) {
  cpblas_pscopy(n, X, ix, jx, descX, incX, Y, iy, jy, descY, incY);
}

inline void pcopy(int n, const double* X, int ix, int jx, const int* descX, int incX,
                  double* Y, int iy, int jy, const int* descY, int incY) {
  cpblas_pdcopy(n, X, ix, jx, descX, incX, Y, iy, jy, descY, incY);
}

inline void pcopy(int n, const std::complex<float>* X, int ix, int jx, const int* descX, int incX,
                  std::complex<float>* Y, int iy, int jy, const int* descY, int incY) {
  cpblas_pccopy(n, complex_cast(X), ix, jx, descX, incX, complex_cast(Y), iy, jy, descY, incY);
}

inline void pcopy(int n, const std::complex<double>* X, int ix, int jx, const int* descX, int incX,
                  std::complex<double>* Y, int iy, int jy, const int* descY, int incY) {
  cpblas_pzcopy(n, complex_cast(X), ix, jx, descX, incX, complex_cast(Y), iy, jy, descY, incY);
}

} // end namespace pblas
} // end namespace rokko
