/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/cpblas.h>

#define CPBLAS_PCOPY_IMPL(NAMES, NAMEL, TYPE) \
void cpblas_ ## NAMES (int N, const TYPE * X, int ix, int jx, const int* DESCX, int INCX, TYPE * Y, int iy, int jy, const int* DESCY, int INCY) { \
  int ix_f = ix + 1; \
  int jx_f = jx + 1; \
  int iy_f = iy + 1; \
  int jy_f = jy + 1; \
  ROKKO_GLOBAL(NAMES, NAMEL) (&N, X, &ix_f, &jx_f, DESCX, &INCX, Y, &iy_f, &jy_f, DESCY, &INCY); \
}

CPBLAS_PCOPY_IMPL(pscopy, PSCOPY, float);
CPBLAS_PCOPY_IMPL(pdcopy, PDCOPY, double);
CPBLAS_PCOPY_IMPL(pccopy, PCCOPY, lapack_complex_float);
CPBLAS_PCOPY_IMPL(pzcopy, PZCOPY, lapack_complex_double);

#undef CPBLAS_PCOPY_IMPL
