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

#define PBLAS_PDOT_IMPL(NAMES, NAMEL, TYPE) \
TYPE PBLAS_ ## NAMES (int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY) { \
  TYPE DOT; ROKKO_GLOBAL(NAMES, NAMEL) (&N, &DOT, X, &IX, &JX, DESCX, &INCX, Y, &IY, &JY, DESCY, &INCY); return DOT; \
} \
void PBLAS_ ## NAMES ## _sub (int N, TYPE * DOT, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY) { \
  ROKKO_GLOBAL(NAMES, NAMEL) (&N, DOT, X, &IX, &JX, DESCX, &INCX, Y, &IY, &JY, DESCY, &INCY); \
}

PBLAS_PDOT_IMPL(psdot, PSDOT, float);
PBLAS_PDOT_IMPL(pddot, PDDOT, double);

#undef PBLAS_PDOT_IMPL
