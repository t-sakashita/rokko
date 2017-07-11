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

#ifndef ROKKO_PBLAS_H
#define ROKKO_PBLAS_H

#include <lapacke.h>
#include <rokko/mangling.h>

/* PvCOPY */

#define PBLAS_PCOPY_DECL(NAMES, NAMEL, TYPE) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const int* N, const TYPE * X, const int* IX, const int* JX, int* DESCX, const int* INCX, TYPE * Y, const int* IY, int* JY, int *DESCY, const int* INCY); \
void PBLAS_ ## NAMES (int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, TYPE * Y, int IY, int JY, int *DESCY, int INCY);

#ifdef __cplusplus
extern "C" {
#endif

PBLAS_PCOPY_DECL(pscopy, PSCOPY, float);
PBLAS_PCOPY_DECL(pdcopy, PDCOPY, double);
PBLAS_PCOPY_DECL(pccopy, PCCOPY, lapack_complex_float);
PBLAS_PCOPY_DECL(pzcopy, PZCOPY, lapack_complex_double);

#ifdef __cplusplus
}
#endif

#undef PBLAS_PCOPY_DECL

/* PvDOT, PvDOTU, PvDOTC */

#define PBLAS_PDOT_DECL(NAMES, NAMEL, TYPE) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const int* N, TYPE * DOT, const TYPE * X, const int* IX, const int* JX, int* DESCX, const int* INCX, const TYPE * Y, const int* IY, const int* JY, int* DESCY, const int* INCY); \
TYPE PBLAS_ ## NAMES (int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY); \
void PBLAS_ ## NAMES ## _sub (int N, TYPE * DOT, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY);

#ifdef __cplusplus
extern "C" {
#endif

PBLAS_PDOT_DECL(psdot, PSDOT, float);
PBLAS_PDOT_DECL(pddot, PDDOT, double);
PBLAS_PDOT_DECL(pcdotu, PCDOTU, lapack_complex_float);
PBLAS_PDOT_DECL(pzdotu, PZDOTU, lapack_complex_double);
PBLAS_PDOT_DECL(pcdotc, PCDOTC, lapack_complex_float);
PBLAS_PDOT_DECL(pzdotc, PZDOTC, lapack_complex_double);

#ifdef __cplusplus
}
#endif

#undef PBLAS_PDOT_DECL

/* PvGEMV */

#define PBLAS_PGEMV_DECL(NAMES, NAMEL, TYPE) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const char* TRANS, const int* M, const int* N, const TYPE * ALPHA, const TYPE * A, const int* IA, const int* JA, int* DESCA, const TYPE * X, const int* IX, const int* JX, int* DESCX, const int* INCX, const TYPE * BETA, TYPE * Y, const int* IY, const int* JY, int* DESCY, const int* INCY); \
void PBLAS_ ## NAMES (char TRANS, int M, int N, TYPE ALPHA, const TYPE * A, int IA, int JA, int* DESCA, const TYPE * X, int IX, int JX, int* DESCX, int INCX, TYPE BETA, TYPE * Y, int IY, int JY, int* DESCY, int INCY);

#ifdef __cplusplus
extern "C" {
#endif

PBLAS_PGEMV_DECL(psgemv, PSGEMV, float);
PBLAS_PGEMV_DECL(pdgemv, PDGEMV, double);
PBLAS_PGEMV_DECL(pcgemv, PCGEMV, lapack_complex_float);
PBLAS_PGEMV_DECL(pzgemv, PZGEMV, lapack_complex_double);

#ifdef __cplusplus
}
#endif

#undef PBLAS_PGEMV_DECL

/* PvGEMM */

#define PBLAS_PGEMM_DECL(NAMES, NAMEL, TYPE) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const char* TRANSA, const char* TRANSB, const int* M, const int* N, const int* K, const TYPE * ALPHA, const TYPE * A, const int* IA, const int* JA, int* DESCA, const TYPE * B, const int* IB, const int* JB, int* DESCB, const TYPE * BETA, TYPE * C, const int* IC, const int* JC, int* DESCC); \
void PBLAS_ ## NAMES (char TRANSA, char TRANSB, int M, int N, int K, TYPE ALPHA, const TYPE * A, int IA, int JA, int* DESCA, const TYPE * B, int IB, int JB, int* DESCB, TYPE BETA, TYPE * C, int IC, int JC, int* DESCC);

#ifdef __cplusplus
extern "C" {
#endif

PBLAS_PGEMM_DECL(psgemm, PSGEMM, float);
PBLAS_PGEMM_DECL(pdgemm, PDGEMM, double);
PBLAS_PGEMM_DECL(pcgemm, PCGEMM, lapack_complex_float);
PBLAS_PGEMM_DECL(pzgemm, PZGEMM, lapack_complex_double);

#ifdef __cplusplus
}
#endif

#undef PBLAS_PGEMM_DECL

#endif // ROKKO_PBLAS_H
