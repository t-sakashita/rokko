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

#ifndef ROKKO_CPBLAS_H
#define ROKKO_CPBLAS_H

#include <rokko/pblas/pblas_interface.h>

/* PvCOPY */

#define CPBLAS_PCOPY_DECL(NAMES, NAMEL, TYPE) \
void cpblas_ ## NAMES (int N, const TYPE * X, int IX, int JX, const int* DESCX, int INCX, TYPE * Y, int IY, int JY, const int* DESCY, int INCY);

#ifdef __cplusplus
extern "C" {
#endif

CPBLAS_PCOPY_DECL(pscopy, PSCOPY, float);
CPBLAS_PCOPY_DECL(pdcopy, PDCOPY, double);
CPBLAS_PCOPY_DECL(pccopy, PCCOPY, lapack_complex_float);
CPBLAS_PCOPY_DECL(pzcopy, PZCOPY, lapack_complex_double);

#ifdef __cplusplus
}
#endif

#undef CPBLAS_PCOPY_DECL

/* PvDOT, PvDOTU, PvDOTC */

#define CPBLAS_PDOT_DECL(NAMES, NAMEL, TYPE) \
TYPE cpblas_ ## NAMES (int N, const TYPE * X, int IX, int JX, const int* DESCX, int INCX, const TYPE * Y, int IY, int JY, const int* DESCY, int INCY); \
void cpblas_ ## NAMES ## _sub (int N, TYPE * DOT, const TYPE * X, int IX, int JX, const int* DESCX, int INCX, const TYPE * Y, int IY, int JY, const int* DESCY, int INCY);

#ifdef __cplusplus
extern "C" {
#endif

CPBLAS_PDOT_DECL(psdot, PSDOT, float);
CPBLAS_PDOT_DECL(pddot, PDDOT, double);
CPBLAS_PDOT_DECL(pcdotu, PCDOTU, lapack_complex_float);
CPBLAS_PDOT_DECL(pzdotu, PZDOTU, lapack_complex_double);
CPBLAS_PDOT_DECL(pcdotc, PCDOTC, lapack_complex_float);
CPBLAS_PDOT_DECL(pzdotc, PZDOTC, lapack_complex_double);

#ifdef __cplusplus
}
#endif

#undef CPBLAS_PDOT_DECL

/* PvGEMV */

#define CPBLAS_PGEMV_DECL(NAMES, NAMEL, TYPE) \
void cpblas_ ## NAMES (char TRANS, int M, int N, TYPE ALPHA, const TYPE * A, int IA, int JA, const int* DESCA, const TYPE * X, int IX, int JX, const int* DESCX, int INCX, TYPE BETA, TYPE * Y, int IY, int JY, const int* DESCY, int INCY);

#ifdef __cplusplus
extern "C" {
#endif

CPBLAS_PGEMV_DECL(psgemv, PSGEMV, float);
CPBLAS_PGEMV_DECL(pdgemv, PDGEMV, double);
CPBLAS_PGEMV_DECL(pcgemv, PCGEMV, lapack_complex_float);
CPBLAS_PGEMV_DECL(pzgemv, PZGEMV, lapack_complex_double);

#ifdef __cplusplus
}
#endif

#undef CPBLAS_PGEMV_DECL

/* PvGEMM */

#define CPBLAS_PGEMM_DECL(NAMES, NAMEL, TYPE) \
void cpblas_ ## NAMES (char TRANSA, char TRANSB, int M, int N, int K, TYPE ALPHA, const TYPE * A, int IA, int JA, const int* DESCA, const TYPE * B, int IB, int JB, const int* DESCB, TYPE BETA, TYPE * C, int IC, int JC, const int* DESCC);

#ifdef __cplusplus
extern "C" {
#endif

CPBLAS_PGEMM_DECL(psgemm, PSGEMM, float);
CPBLAS_PGEMM_DECL(pdgemm, PDGEMM, double);
CPBLAS_PGEMM_DECL(pcgemm, PCGEMM, lapack_complex_float);
CPBLAS_PGEMM_DECL(pzgemm, PZGEMM, lapack_complex_double);

#ifdef __cplusplus
}
#endif

#undef CPBLAS_PGEMM_DECL

#endif // ROKKO_CPBLAS_H
