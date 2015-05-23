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

#define PBLAS_PDOT_IMPL(NAMES, NAMEL, TYPEC, TYPEX) \
extern "C" { \
void PBLASE_ ## NAMES (int N, TYPEC * DOT, const TYPEC * X, int IX, int JX, int* DESCX, int INCX, const TYPEC * Y, int IY, int JY, int* DESCY, int INCY) { \
  LAPACK_GLOBAL(NAMES, NAMEL) (&N, DOT, X, &IX, &JX, DESCX, &INCX, Y, &IY, &JY, DESCY, &INCY); } \
} \
void PBLASE_pdot(int N, TYPEC * DOT, const TYPEX * X, int IX, int JX, int* DESCX, int INCX, const TYPEX * Y, int IY, int JY, int* DESCY, int INCY) { \
  PBLASE_ ## NAMES (N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY); }

PBLAS_PDOT_IMPL(psdot, PSDOT, float, float);
PBLAS_PDOT_IMPL(pddot, PDDOT, double, double);

#undef PBLAS_PDOT_IMPL
