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

#include <rokko/mangling.h>

/* PvCOPY */

#ifdef __cplusplus

#define PBLAS_PCOPY_DECL(NAMES, NAMEL, TYPEC, TYPEX) \
extern "C" { \
void LAPACK_GLOBAL(NAMES, NAMEL) (const int* N, const TYPEC * X, const int* IX, const int* JX, int* DESCX, const int* INCX, TYPEC * Y, const int* IY, const int* JY, int *DESCY, const int* INCY); \
void PBLASE_ ## NAMES (int N, const TYPEC * X, int IX, int JX, int* DESCX, int INCX, TYPEC * Y, int IY, int JY, int *DESCY, int INCY); \
} \
void PBLASE_pcopy(int N, const TYPEX * X, int IX, int JX, int* DESCX, int INCX, TYPEX * Y, int IY, int JY, int *DESCY, int INCY);

#else

#define PBLAS_PCOPY_DECL(NAMES, NAMEL, TYPEC, TYPEX) \
void LAPACK_GLOBAL(NAMES, NAMEL) (const int* N, const TYPEC * X, const int* IX, const int* JX, int* DESCX, const int* INCX, TYPEC * Y, const int* IY, const int* JY, int *DESCY, const int* INCY); \
void PBLASE_ ## NAMES (int N, const TYPEC * X, int IX, int JX, int* DESCX, int INCX, TYPEC * Y, int IY, int JY, int *DESCY, int INCY); \

#endif

PBLAS_PCOPY_DECL(pscopy, PSCOPY, float, float);
PBLAS_PCOPY_DECL(pdcopy, PDCOPY, double, double);
PBLAS_PCOPY_DECL(pccopy, PCCOPY, lapack_complex_float, std::complex<float>);
PBLAS_PCOPY_DECL(pzcopy, PZCOPY, lapack_complex_double, std::complex<double>);

#undef PBLAS_PCOPY_DECL

/* PvDOT */

#ifdef __cplusplus

#define PBLAS_PDOT_DECL(NAMES, NAMEL, TYPEC, TYPEX) \
extern "C" { \
void LAPACK_GLOBAL(NAMES, NAMEL) (const int* N, TYPEC * DOT, const TYPEC * X, const int* IX, const int* JX, int* DESCX, const int* INCX, const TYPEC * Y, const int* IY, const int* JY, int* DESCY, const int* INCY); \
void PBLASE_ ## NAMES (int N, TYPEC * DOT, const TYPEC * X, int IX, int JX, int* DESCX, int INCX, const TYPEC * Y, int IY, int JY, int* DESCY, int INCY); \
} \
void PBLASE_pdot(int N, TYPEC * DOT, const TYPEX * X, int IX, int JX, int* DESCX, int INCX, const TYPEX * Y, int IY, int JY, int* DESCY, int INCY);

#else

#define PBLAS_PDOT_DECL(NAMES, NAMEL, TYPEC, TYPEX) \
void LAPACK_GLOBAL(NAMES, NAMEL) (const int* N, TYPEC * DOT, const TYPEC * X, const int* IX, const int* JX, int* DESCX, const int* INCX, const TYPEC * Y, const int* IY, const int* JY, int* DESCY, const int* INCY); \
void PBLASE_ ## NAMES (int N, TYPEC * DOT, const TYPEC * X, int IX, int JX, int* DESCX, int INCX, const TYPEC * Y, int IY, int JY, int* DESCY, int INCY); \

#endif

PBLAS_PDOT_DECL(psdot, PSDOT, float, float);
PBLAS_PDOT_DECL(pddot, PDDOT, double, double);

#undef PBLAS_PVDOT

/* PvGEMV */

#ifdef __cplusplus

#define PBLAS_PGEMV_DECL(NAMES, NAMEL, TYPEC, TYPEX) \
extern "C" { \
void LAPACK_GLOBAL(NAMES, NAMEL) (const char* TRANS, const int* M, const int* N, const TYPEC * ALPHA, const TYPEC * A, const int* IA, const int* JA, int* DESCA, const TYPEC * X, const int* IX, const int* JX, int* DESCX, const int* INCX, const TYPEC * BETA, TYPEC * Y, const int* IY, const int* JY, int* DESCY, const int* INCY); \
void PBLASE_ ## NAMES (char TRANS, int M, int N, TYPEC ALPHA, const TYPEC * A, int IA, int JA, int* DESCA, const TYPEC * X, int IX, int JX, int* DESCX, int INCX, TYPEC BETA, TYPEC * Y, int IY, int JY, int* DESCY, int INCY); \
} \
void PBLASE_pgemv(char TRANS, int M, int N, TYPEX ALPHA, const TYPEX * A, int IA, int JA, int* DESCA, const TYPEX * X, int IX, int JX, int* DESCX, int INCX, TYPEX BETA, TYPEX * Y, int IY, int JY, int* DESCY, int INCY);

#else

#define PBLAS_PGEMV_DECL(NAMES, NAMEL, TYPEC, TYPEX) \
void LAPACK_GLOBAL(NAMES, NAMEL) (const char* TRANS, const int* M, const int* N, const TYPEC * ALPHA, const TYPEC * A, const int* IA, const int* JA, int* DESCA, const TYPEC * X, const int* IX, const int* JX, int* DESCX, const int* INCX, const TYPEC * BETA, TYPEC * Y, const int* IY, const int* JY, int* DESCY, const int* INCY); \
void PBLASE_ ## NAMES (char TRANS, int M, int N, TYPEC ALPHA, const TYPEC * A, int IA, int JA, int* DESCA, const TYPEC * X, int IX, int JX, int* DESCX, int INCX, TYPEC BETA, TYPEC * Y, int IY, int JY, int* DESCY, int INCY);

#endif

PBLAS_PGEMV_DECL(psgemv, PSGEMV, float, float);
PBLAS_PGEMV_DECL(pdgemv, PDGEMV, double, double);
PBLAS_PGEMV_DECL(pcgemv, PCGEMV, lapack_complex_float, std::complex<float>);
PBLAS_PGEMV_DECL(pzgemv, PZGEMV, lapack_complex_double, std::complex<double>);

#undef PBLAS_PVGEMV

/* PvGEMM */

#ifdef __cplusplus

#define PBLAS_PGEMM_DECL(NAMES, NAMEL, TYPEC, TYPEX) \
extern "C" { \
void LAPACK_GLOBAL(NAMES, NAMEL) (const char* TRANSA, const char* TRANSB, const int* M, const int* N, const int* K, const TYPEC * ALPHA, const TYPEC * A, const int* IA, const int* JA, int* DESCA, const TYPEC * B, const int* IB, const int* JB, int* DESCB, const TYPEC * BETA, TYPEC * C, const int* IC, const int* JC, int* DESCC); \
void PBLASE_ ## NAMES (char TRANSA, char TRANSB, int M, int N, int K, TYPEC ALPHA, const TYPEC * A, int IA, int JA, int* DESCA, const TYPEC * B, int IB, int JB, int* DESCB, TYPEC BETA, TYPEC * C, int IC, int JC, int* DESCC); \
} \
void PBLASE_pgemm(char TRANSA, char TRANSB, int M, int N, int K, TYPEX ALPHA, const TYPEX * A, int IA, int JA, int* DESCA, const TYPEX * B, int IB, int JB, int* DESCB, TYPEX BETA, TYPEX * C, int IC, int JC, int* DESCC);

#else

#define PBLAS_PGEMM_DECL(NAMES, NAMEL, TYPEC, TYPEX) \
void LAPACK_GLOBAL(NAMES, NAMEL) (const char* TRANSA, const char* TRANSB, const int* M, const int* N, const int* K, const TYPEC * ALPHA, const TYPEC * A, const int* IA, const int* JA, int* DESCA, const TYPEC * B, const int* IB, const int* JB, int* DESCB, const TYPEC * BETA, TYPEC * C, const int* IC, const int* JC, int* DESCC); \
void PBLASE_ ## NAMES (char TRANSA, char TRANSB, int M, int N, int K, TYPEC ALPHA, const TYPEC * A, int IA, int JA, int* DESCA, const TYPEC * B, int IB, int JB, int* DESCB, TYPEC BETA, TYPEC * C, int IC, int JC, int* DESCC);

#endif

PBLAS_PGEMM_DECL(psgemm, PSGEMM, float, float);
PBLAS_PGEMM_DECL(pdgemm, PDGEMM, double, double);
PBLAS_PGEMM_DECL(pcgemm, PCGEMM, lapack_complex_float, std::complex<float>);
PBLAS_PGEMM_DECL(pzgemm, PZGEMM, lapack_complex_double, std::complex<double>);

#undef PBLAS_PGEMM_DECL

#endif // ROKKO_PBLAS_H
