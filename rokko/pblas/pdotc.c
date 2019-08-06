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

#define CPBLAS_PDOTC_IMPL(NAMES, NAMEL, TYPE) \
TYPE cpblas_ ## NAMES (int n, const TYPE * x, int ix, int jx, const int* descx, int incx, const TYPE * y, int iy, int jy, const int* descy, int incy) { \
  int ix_f = ix + 1; \
  int jx_f = jx + 1; \
  int iy_f = iy + 1; \
  int jy_f = jy + 1; \
  TYPE dotc; \
  ROKKO_GLOBAL(NAMES, NAMEL) (&n, &dotc, x, &ix_f, &jx_f, descx, &incx, y, &iy_f, &jy_f, descy, &incy); \
  return dotc; \
} \
void cpblas_ ## NAMES ## _sub (int n, TYPE * dotc, const TYPE * x, int ix, int jx, const int* descx, int incx, const TYPE * y, int iy, int jy, const int* descy, int incy) { \
  int ix_f = ix + 1; \
  int jx_f = jx + 1; \
  int iy_f = iy + 1; \
  int jy_f = jy + 1; \
  ROKKO_GLOBAL(NAMES, NAMEL) (&n, dotc, x, &ix_f, &jx_f, descx, &incx, y, &iy_f, &jy_f, descy, &incy); \
}

CPBLAS_PDOTC_IMPL(pcdotc, PCDOTC, lapack_complex_float);
CPBLAS_PDOTC_IMPL(pzdotc, PZDOTC, lapack_complex_double);

#undef CPBLAS_PDOTC_IMPL
