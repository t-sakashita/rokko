/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/cscalapack.h>
#include <rokko/scalapack/scalapack_interface.h>

#define CSCALAPACK_PSYEVD_WORK_IMPL(NAMES, NAMEL, TYPE) \
int cscalapack_ ## NAMES ## _work(char jobz, char uplo, int n, \
                                  TYPE* A, int ia, int ja, const int* descA, \
                                  TYPE* w, TYPE* Z, int iz, int jz, const int* descZ, \
                                  TYPE* work, int lwork, int* iwork, int liwork) { \
  int ia_f = ia + 1; \
  int ja_f = ja + 1; \
  int iz_f = iz + 1; \
  int jz_f = jz + 1; \
  int info; \
  ROKKO_GLOBAL(NAMES, NAMEL) (&jobz, &uplo, &n, A, &ia_f, &ja_f, descA, w, Z, &iz_f, &jz_f, descZ, \
                              work, &lwork, iwork, &liwork, &info); \
  return info; \
}

CSCALAPACK_PSYEVD_WORK_IMPL(pssyevd, PSSYEVD, float)
CSCALAPACK_PSYEVD_WORK_IMPL(pdsyevd, PDSYEVD, double)
CSCALAPACK_PSYEVD_WORK_IMPL(pcheevd, PCHEEVD, lapack_complex_float)
CSCALAPACK_PSYEVD_WORK_IMPL(pzheevd, PZHEEVD, lapack_complex_double)

#undef CSCALAPACK_PSYEVD_WORK_IMPL
