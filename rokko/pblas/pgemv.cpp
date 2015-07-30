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

#include <rokko/pblas.h>

#define PBLAS_PGEMV_IMPL(NAMES, NAMEL, TYPEC, TYPEX) \
extern "C" { \
void PBLASE_ ## NAMES (char TRANS, int M, int N, TYPEC ALPHA, const TYPEC * A, int IA, int JA, int* DESCA, const TYPEC * X, int IX, int JX, int* DESCX, int INCX, TYPEC BETA, TYPEC * Y, int IY, int JY, int* DESCY, int INCY) { \
  LAPACK_GLOBAL(NAMES, NAMEL) (&TRANS, &M, &N, &ALPHA, A, &IA, &JA, DESCA, X, &IX, &JX, DESCX, &INCX, &BETA, Y, &IY, &JY, DESCY, &INCY); } \
} \
void PBLASE_pgemv(char TRANS, int M, int N, TYPEX ALPHA, const TYPEX * A, int IA, int JA, int* DESCA, const TYPEX * X, int IX, int JX, int* DESCX, int INCX, TYPEX BETA, TYPEX * Y, int IY, int JY, int* DESCY, int INCY) { \
  PBLASE_ ## NAMES (TRANS, M, N, ALPHA, A, IA, JA, DESCA, X, IX, JX, DESCX, INCX, BETA, Y, IY, JY, DESCY, INCY); }

PBLAS_PGEMV_IMPL(psgemv, PSGEMV, float, float);
PBLAS_PGEMV_IMPL(pdgemv, PDGEMV, double, double);
PBLAS_PGEMV_IMPL(pcgemv, PCGEMV, lapack_complex_float, std::complex<float>);
PBLAS_PGEMV_IMPL(pzgemv, PZGEMV, lapack_complex_double, std::complex<double>);

#undef PBLAS_PGEMV_IMPL
