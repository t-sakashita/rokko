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

#define CSCALAPACK_PSYEV_WORK_IMPL(NAMES, NAMEL, TYPE, TYPE_REAL) \
int cscalapack_ ## NAMES ## _work(char jobz, char uplo, int n, \
                                  TYPE* A, int ia, int ja, const int* descA, \
                                  TYPE_REAL* w, TYPE* Z, int iz, int jz, const int* descZ, \
                                  TYPE* work, int lwork) { \
  int ia_f = ia + 1; \
  int ja_f = ja + 1; \
  int iz_f = iz + 1; \
  int jz_f = jz + 1; \
  int info; \
  ROKKO_GLOBAL(NAMES, NAMEL) (&jobz, &uplo, &n, A, &ia_f, &ja_f, descA, w, Z, &iz_f, &jz_f, descZ, \
                              work, &lwork, &info); \
  return info; \
}

CSCALAPACK_PSYEV_WORK_IMPL(pssyev, PSSYEV, float, float)
CSCALAPACK_PSYEV_WORK_IMPL(pdsyev, PDSYEV, double, double)
CSCALAPACK_PSYEV_WORK_IMPL(pcheev, PCHEEV, lapack_complex_float, float)
CSCALAPACK_PSYEV_WORK_IMPL(pzheev, PZHEEV, lapack_complex_double, double)

#undef CSCALAPACK_PSYEV_WORK_IMPL
