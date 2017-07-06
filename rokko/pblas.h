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
#ifdef __cplusplus
# include <boost/numeric/ublas/traits.hpp>
#endif

/* PvCOPY */

#ifdef __cplusplus

#define PBLAS_PCOPY_DECL(NAMES, NAMEL, TYPE) \
extern "C" { \
void ROKKO_GLOBAL(NAMES, NAMEL) (const int* N, const TYPE * X, const int* IX, const int* JX, int* DESCX, const int* INCX, TYPE * Y, const int* IY, const int* JY, int *DESCY, const int* INCY); \
void PBLASE_ ## NAMES (int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, TYPE * Y, int IY, int JY, int *DESCY, int INCY); \
} \
void PBLASE_pcopy(int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, TYPE * Y, int IY, int JY, int *DESCY, int INCY);

#else

#define PBLAS_PCOPY_DECL(NAMES, NAMEL, TYPE) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const int* N, const TYPE * X, const int* IX, const int* JX, int* DESCX, const int* INCX, TYPE * Y, const int* IY, const int* JY, int *DESCY, const int* INCY); \
void PBLASE_ ## NAMES (int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, TYPE * Y, int IY, int JY, int *DESCY, int INCY); \

#endif

PBLAS_PCOPY_DECL(pscopy, PSCOPY, float);
PBLAS_PCOPY_DECL(pdcopy, PDCOPY, double);
PBLAS_PCOPY_DECL(pccopy, PCCOPY, lapack_complex_float);
PBLAS_PCOPY_DECL(pzcopy, PZCOPY, lapack_complex_double);

#undef PBLAS_PCOPY_DECL

/* PvDOT */

#ifdef __cplusplus

#define PBLAS_PDOT_DECL(NAMES, NAMEL, TYPE) \
extern "C" { \
void ROKKO_GLOBAL(NAMES, NAMEL) (const int* N, TYPE * DOT, const TYPE * X, const int* IX, const int* JX, int* DESCX, const int* INCX, const TYPE * Y, const int* IY, const int* JY, int* DESCY, const int* INCY); \
TYPE PBLASE_ ## NAMES (int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY); \
} \
TYPE PBLASE_pdot(int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY); \
TYPE PBLASE_pdotc(int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY);

#else

#define PBLAS_PDOT_DECL(NAMES, NAMEL, TYPE) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const int* N, TYPE * DOT, const TYPE * X, const int* IX, const int* JX, int* DESCX, const int* INCX, const TYPE * Y, const int* IY, const int* JY, int* DESCY, const int* INCY); \
TYPE PBLASE_ ## NAMES (int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY); \

#endif

PBLAS_PDOT_DECL(psdot, PSDOT, float);
PBLAS_PDOT_DECL(pddot, PDDOT, double);

/* PvDOTC */

#ifdef __cplusplus

#undef PBLAS_PDOT_DECL

#define PBLAS_PDOTC_DECL(NAMES, NAMEL, TYPE) \
extern "C" { \
void ROKKO_GLOBAL(NAMES, NAMEL) (const int* N, boost::numeric::ublas::type_traits<TYPE>::real_type * DOTC, const TYPE * X, const int* IX, const int* JX, int* DESCX, const int* INCX, const TYPE * Y, const int* IY, const int* JY, int* DESCY, const int* INCY); \
boost::numeric::ublas::type_traits<TYPE>::real_type PBLASE_ ## NAMES (int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY); \
} \
boost::numeric::ublas::type_traits<TYPE>::real_type PBLASE_pdotc(int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY);

#else

void ROKKO_GLOBAL(NAMES, NAMEL) (const int* N, boost::numeric::ublas::type_traits<TYPE>::real_type * DOTC, const TYPE * X, const int* IX, const int* JX, int* DESCX, const int* INCX, const TYPE * Y, const int* IY, const int* JY, int* DESCY, const int* INCY); \
boost::numeric::ublas::type_traits<TYPE>::real_type PBLASE_ ## NAMES (int N, const TYPE * X, int IX, int JX, int* DESCX, int INCX, const TYPE * Y, int IY, int JY, int* DESCY, int INCY); \

#endif

PBLAS_PDOTC_DECL(pcdotc, PCDOTC, lapack_complex_float);
PBLAS_PDOTC_DECL(pzdotc, PZDOTC, lapack_complex_double);

#undef PBLAS_PDOTC_DECL

/* PvGEMV */

#ifdef __cplusplus

#define PBLAS_PGEMV_DECL(NAMES, NAMEL, TYPE) \
extern "C" { \
void ROKKO_GLOBAL(NAMES, NAMEL) (const char* TRANS, const int* M, const int* N, const TYPE * ALPHA, const TYPE * A, const int* IA, const int* JA, int* DESCA, const TYPE * X, const int* IX, const int* JX, int* DESCX, const int* INCX, const TYPE * BETA, TYPE * Y, const int* IY, const int* JY, int* DESCY, const int* INCY); \
void PBLASE_ ## NAMES (char TRANS, int M, int N, TYPE ALPHA, const TYPE * A, int IA, int JA, int* DESCA, const TYPE * X, int IX, int JX, int* DESCX, int INCX, TYPE BETA, TYPE * Y, int IY, int JY, int* DESCY, int INCY); \
} \
void PBLASE_pgemv(char TRANS, int M, int N, TYPE ALPHA, const TYPE * A, int IA, int JA, int* DESCA, const TYPE * X, int IX, int JX, int* DESCX, int INCX, TYPE BETA, TYPE * Y, int IY, int JY, int* DESCY, int INCY);

#else

#define PBLAS_PGEMV_DECL(NAMES, NAMEL, TYPE) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const char* TRANS, const int* M, const int* N, const TYPE * ALPHA, const TYPE * A, const int* IA, const int* JA, int* DESCA, const TYPE * X, const int* IX, const int* JX, int* DESCX, const int* INCX, const TYPE * BETA, TYPE * Y, const int* IY, const int* JY, int* DESCY, const int* INCY); \
void PBLASE_ ## NAMES (char TRANS, int M, int N, TYPE ALPHA, const TYPE * A, int IA, int JA, int* DESCA, const TYPE * X, int IX, int JX, int* DESCX, int INCX, TYPE BETA, TYPE * Y, int IY, int JY, int* DESCY, int INCY);

#endif

PBLAS_PGEMV_DECL(psgemv, PSGEMV, float);
PBLAS_PGEMV_DECL(pdgemv, PDGEMV, double);
PBLAS_PGEMV_DECL(pcgemv, PCGEMV, lapack_complex_float);
PBLAS_PGEMV_DECL(pzgemv, PZGEMV, lapack_complex_double);

#undef PBLAS_PGEMV_DECL

/* PvGEMM */

#ifdef __cplusplus

#define PBLAS_PGEMM_DECL(NAMES, NAMEL, TYPE) \
extern "C" { \
void ROKKO_GLOBAL(NAMES, NAMEL) (const char* TRANSA, const char* TRANSB, const int* M, const int* N, const int* K, const TYPE * ALPHA, const TYPE * A, const int* IA, const int* JA, int* DESCA, const TYPE * B, const int* IB, const int* JB, int* DESCB, const TYPE * BETA, TYPE * C, const int* IC, const int* JC, int* DESCC); \
void PBLASE_ ## NAMES (char TRANSA, char TRANSB, int M, int N, int K, TYPE ALPHA, const TYPE * A, int IA, int JA, int* DESCA, const TYPE * B, int IB, int JB, int* DESCB, TYPE BETA, TYPE * C, int IC, int JC, int* DESCC); \
} \
void PBLASE_pgemm(char TRANSA, char TRANSB, int M, int N, int K, TYPE ALPHA, const TYPE * A, int IA, int JA, int* DESCA, const TYPE * B, int IB, int JB, int* DESCB, TYPE BETA, TYPE * C, int IC, int JC, int* DESCC);

#else

#define PBLAS_PGEMM_DECL(NAMES, NAMEL, TYPE) \
void ROKKO_GLOBAL(NAMES, NAMEL) (const char* TRANSA, const char* TRANSB, const int* M, const int* N, const int* K, const TYPE * ALPHA, const TYPE * A, const int* IA, const int* JA, int* DESCA, const TYPE * B, const int* IB, const int* JB, int* DESCB, const TYPE * BETA, TYPE * C, const int* IC, const int* JC, int* DESCC); \
void PBLASE_ ## NAMES (char TRANSA, char TRANSB, int M, int N, int K, TYPE ALPHA, const TYPE * A, int IA, int JA, int* DESCA, const TYPE * B, int IB, int JB, int* DESCB, TYPE BETA, TYPE * C, int IC, int JC, int* DESCC);

#endif

PBLAS_PGEMM_DECL(psgemm, PSGEMM, float);
PBLAS_PGEMM_DECL(pdgemm, PDGEMM, double);
PBLAS_PGEMM_DECL(pcgemm, PCGEMM, lapack_complex_float);
PBLAS_PGEMM_DECL(pzgemm, PZGEMM, lapack_complex_double);

#undef PBLAS_PGEMM_DECL

#endif // ROKKO_PBLAS_H
