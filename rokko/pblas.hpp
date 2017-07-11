/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2017 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_PBLAS_HPP
#define ROKKO_PBLAS_HPP

#include <rokko/pblas.h>
#include <rokko/lapack/complex_cast.hpp>

using rokko::lapack::complex_cast;

namespace rokko {
namespace pblas {

// pcopy

#define PBLAS_PCOPY_IMPL(NAMES, TYPE) \
inline void pcopy(int N, const TYPE* X, int IX, int JX, int* DESCX, int INCX, TYPE* Y, int IY, int JY, int *DESCY, int INCY) { \
  PBLAS_ ## NAMES (N, complex_cast(X), IX, JX, DESCX, INCX, complex_cast(Y), IY, JY, DESCY, INCY); \
}

PBLAS_PCOPY_IMPL(pscopy, float);
PBLAS_PCOPY_IMPL(pdcopy, double);
PBLAS_PCOPY_IMPL(pccopy, std::complex<float>);
PBLAS_PCOPY_IMPL(pzcopy, std::complex<double>);

#undef PBLAS_PCOPY_IMPL

// pdot, pdotu, pdotc

#define PBLAS_PDOT_IMPL(NAMEX, NAMES, TYPE) \
inline TYPE NAMEX (int N, const TYPE* X, int IX, int JX, int* DESCX, int INCX, const TYPE* Y, int IY, int JY, int* DESCY, int INCY) { \
  TYPE DOT; \
  PBLAS_ ## NAMES ## _sub (N, complex_cast(&DOT), complex_cast(X), IX, JX, DESCX, INCX, complex_cast(Y), IY, JY, DESCY, INCY); \
  return DOT; \
}

PBLAS_PDOT_IMPL(pdot, psdot, float);
PBLAS_PDOT_IMPL(pdot, pddot, double);
PBLAS_PDOT_IMPL(pdot, pcdotc, std::complex<float>);
PBLAS_PDOT_IMPL(pdot, pzdotc, std::complex<double>);
PBLAS_PDOT_IMPL(pdotu, pcdotu, std::complex<float>);
PBLAS_PDOT_IMPL(pdotu, pzdotu, std::complex<double>);
PBLAS_PDOT_IMPL(pdotc, psdot, float);
PBLAS_PDOT_IMPL(pdotc, pddot, double);
PBLAS_PDOT_IMPL(pdotc, pcdotc, std::complex<float>);
PBLAS_PDOT_IMPL(pdotc, pzdotc, std::complex<double>);

#undef PBLAS_PDOT_IMPL

// pgemv

#define PBLAS_PGEMV_IMPL(NAMES, TYPE) \
inline void pgemv(char TRANS, int M, int N, TYPE ALPHA, const TYPE* A, int IA, int JA, int* DESCA, const TYPE* X, int IX, int JX, int* DESCX, int INCX, TYPE BETA, TYPE* Y, int IY, int JY, int* DESCY, int INCY) { \
  PBLAS_ ## NAMES (TRANS, M, N, complex_cast(ALPHA), complex_cast(A), IA, JA, DESCA, complex_cast(X), IX, JX, DESCX, INCX, complex_cast(BETA), complex_cast(Y), IY, JY, DESCY, INCY); \
}

PBLAS_PGEMV_IMPL(psgemv, float);
PBLAS_PGEMV_IMPL(pdgemv, double);
PBLAS_PGEMV_IMPL(pcgemv, std::complex<float>);
PBLAS_PGEMV_IMPL(pzgemv, std::complex<double>);

#undef PBLAS_PGEMV_IMPL

// pgemm

#define PBLAS_PGEMM_IMPL(NAMES, TYPE) \
inline void pgemm(char TRANSA, char TRANSB, int M, int N, int K, TYPE ALPHA, const TYPE* A, int IA, int JA, int* DESCA, const TYPE* B, int IB, int JB, int* DESCB, TYPE BETA, TYPE* C, int IC, int JC, int* DESCC) { \
  PBLAS_ ## NAMES (TRANSA, TRANSB, M, N, K, complex_cast(ALPHA), complex_cast(A), IA, JA, DESCA, complex_cast(B), IB, JB, DESCB, complex_cast(BETA), complex_cast(C), IC, JC, DESCC); \
}

PBLAS_PGEMM_IMPL(psgemm, float);
PBLAS_PGEMM_IMPL(pdgemm, double);
PBLAS_PGEMM_IMPL(pcgemm, std::complex<float>);
PBLAS_PGEMM_IMPL(pzgemm, std::complex<double>);

#undef PBLAS_PGEMM_IMPL

} // namespace pblas
} // namespace rokko

#endif // ROKKO_PBLAS_HPP
