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

#include <rokko/pblas.h>

#define PBLAS_PCOPY_IMPL(NAMES, NAMEL, TYPE) \
void PBLAS_ ## NAMES (int N, const TYPE * X, int IX, int JX, const int* DESCX, int INCX, TYPE * Y, int IY, int JY, const int* DESCY, int INCY) { \
  ROKKO_GLOBAL(NAMES, NAMEL) (&N, X, &IX, &JX, DESCX, &INCX, Y, &IY, &JY, DESCY, &INCY); }

PBLAS_PCOPY_IMPL(pscopy, PSCOPY, float);
PBLAS_PCOPY_IMPL(pdcopy, PDCOPY, double);
PBLAS_PCOPY_IMPL(pccopy, PCCOPY, lapack_complex_float);
PBLAS_PCOPY_IMPL(pzcopy, PZCOPY, lapack_complex_double);

#undef PBLAS_PCOPY_IMPL
