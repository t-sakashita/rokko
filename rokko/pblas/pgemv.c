/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/cpblas.h>

#define CPBLAS_PGEMV_IMPL(NAMES, NAMEL, TYPEC, TYPEX) \
void cpblas_ ## NAMES (char trans, int m, int n, TYPEC alpha, const TYPEC * a, int ia, int ja, const int* desca, const TYPEC * x, int ix, int jx, const int* descx, int incx, TYPEC beta, TYPEC * y, int iy, int jy, const int* descy, int incy) { \
  int ia_f = ia + 1; \
  int ja_f = ja + 1; \
  int ix_f = ix + 1; \
  int jx_f = jx + 1; \
  int iy_f = iy + 1; \
  int jy_f = jy + 1; \
  ROKKO_GLOBAL(NAMES, NAMEL) (&trans, &m, &n, &alpha, a, &ia_f, &ja_f, desca, x, &ix_f, &jx_f, descx, &incx, &beta, y, &iy_f, &jy_f, descy, &incy); \
}

CPBLAS_PGEMV_IMPL(psgemv, PSGEMV, float, float);
CPBLAS_PGEMV_IMPL(pdgemv, PDGEMV, double, double);
CPBLAS_PGEMV_IMPL(pcgemv, PCGEMV, lapack_complex_float, std::complex<float>);
CPBLAS_PGEMV_IMPL(pzgemv, PZGEMV, lapack_complex_double, std::complex<double>);

#undef CPBLAS_PGEMV_IMPL
