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

#define PBLAS_PDOTU_IMPL(NAMES, NAMEL, TYPE) \
TYPE PBLAS_ ## NAMES (int N, const TYPE * X, int IX, int JX, const int* DESCX, int INCX, const TYPE * Y, int IY, int JY, const int* DESCY, int INCY) { \
  TYPE DOTU; ROKKO_GLOBAL(NAMES, NAMEL) (&N, &DOTU, X, &IX, &JX, DESCX, &INCX, Y, &IY, &JY, DESCY, &INCY); return DOTU; \
} \
void PBLAS_ ## NAMES ## _sub (int N, TYPE * DOTU, const TYPE * X, int IX, int JX, const int* DESCX, int INCX, const TYPE * Y, int IY, int JY, const int* DESCY, int INCY) { \
  ROKKO_GLOBAL(NAMES, NAMEL) (&N, DOTU, X, &IX, &JX, DESCX, &INCX, Y, &IY, &JY, DESCY, &INCY); \
}

PBLAS_PDOTU_IMPL(pcdotu, PCDOTU, lapack_complex_float);
PBLAS_PDOTU_IMPL(pzdotu, PZDOTU, lapack_complex_double);

#undef PBLAS_PDOTU_IMPL
