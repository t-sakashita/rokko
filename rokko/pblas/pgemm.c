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

#include <rokko/cpblas.h>
#include <rokko/pblas/pblas_interface.h>

#define CPBLAS_PGEMM_IMPL(NAMES, NAMEL, TYPEC, TYPEX) \
void cpblas_ ## NAMES (char transa, char transb, int m, int n, int k, TYPEC alpha, const TYPEC * a, int ia, int ja, const int* desca, const TYPEC * b, int ib, int jb, const int* descb, TYPEC beta, TYPEC * c, int ic, int jc, const int* descc) { \
  int ia_f = ia + 1; \
  int ja_f = ja + 1; \
  int ib_f = ib + 1; \
  int jb_f = jb + 1; \
  int ic_f = ic + 1; \
  int jc_f = jc + 1; \
  ROKKO_GLOBAL(NAMES, NAMEL) (&transa, &transb, &m, &n, &k, &alpha, a, &ia_f, &ja_f, desca, b, &ib_f, &jb_f, descb, &beta, c, &ic_f, &jc_f, descc); \
}

CPBLAS_PGEMM_IMPL(psgemm, PSGEMM, float, float);
CPBLAS_PGEMM_IMPL(pdgemm, PDGEMM, double, double);
CPBLAS_PGEMM_IMPL(pcgemm, PCGEMM, lapack_complex_float, std::complex<float>);
CPBLAS_PGEMM_IMPL(pzgemm, PZGEMM, lapack_complex_double, std::complex<double>);

#undef CPBLAS_PGEMM_IMPL
