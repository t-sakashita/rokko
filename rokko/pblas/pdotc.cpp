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
TYPE PBLASE_pdotc(int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY) { \
  return PBLASE_ ## NAMES (N, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY); }

PBLAS_PDOT_IMPL(psdot, PSDOT, float);
PBLAS_PDOT_IMPL(pddot, PDDOT, double);

#undef PBLAS_PDOT_IMPL

#define PBLAS_PDOTC_IMPL(NAMES, NAMEL, TYPE) \
extern "C" { \
boost::numeric::ublas::type_traits<TYPE>::real_type PBLASE_ ## NAMES (int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY) { \
  boost::numeric::ublas::type_traits<TYPE>::real_type DOTC; LAPACK_GLOBAL(NAMES, NAMEL) (&N, &DOTC, X, &IX, &JX, DESCX, &INCX, Y, &IY, &JY, DESCY, &INCY); return DOTC; } \
} \
boost::numeric::ublas::type_traits<TYPE>::real_type PBLASE_pdotc(int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY) { \
  return PBLASE_ ## NAMES (N, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY); }
  
PBLAS_PDOTC_IMPL(pcdotc, PCDOTC, lapack_complex_float);
PBLAS_PDOTC_IMPL(pzdotc, PZDOTC, lapack_complex_double);

#undef PBLAS_PDOTC_IMPL
