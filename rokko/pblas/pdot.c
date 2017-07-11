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
extern "C" {                                                          \
TYPE PBLASE_ ## NAMES (int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY) { \
  TYPE DOT; ROKKO_GLOBAL(NAMES, NAMEL) (&N, &DOT, X, &IX, &JX, DESCX, &INCX, Y, &IY, &JY, DESCY, &INCY); return DOT; } \
} \
TYPE PBLASE_pdot(int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY) { \
  return PBLASE_ ## NAMES (N, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY); }

PBLAS_PDOT_IMPL(psdot, PSDOT, float);
PBLAS_PDOT_IMPL(pddot, PDDOT, double);

#undef PBLAS_PDOT_IMPL
