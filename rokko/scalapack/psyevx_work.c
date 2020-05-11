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

#define CSCALAPACK_PSYEVX_WORK_IMPL(NAMES, NAMEL, TYPE, TYPE_REAL) \
int cscalapack_ ## NAMES ## _work(char jobz, char range, char uplo, int n, \
                                  TYPE* A, int ia, int ja, const int* descA, \
                                  TYPE_REAL vl, TYPE_REAL vu, int il, int iu, \
                                  TYPE_REAL abstol, int* m, int* nZ, TYPE_REAL* w, TYPE_REAL orfac, \
                                  TYPE* Z, int iz, int jz, const int* descZ, \
                                  TYPE* work, int lwork, int* iwork, int liwork, \
                                  int* ifail, int* iclustr, TYPE_REAL* gap) { \
  int ia_f = ia + 1; \
  int ja_f = ja + 1; \
  int il_f = il + 1; \
  int iu_f = iu + 1; \
  int iz_f = iz + 1; \
  int jz_f = jz + 1; \
  int info; \
  ROKKO_GLOBAL(NAMES, NAMEL) (&jobz, &range, &uplo, &n, A, &ia_f, &ja_f, descA, &vl, &vu, &il_f, &iu_f, \
                              &abstol, m, nZ, w, &orfac, Z, &iz_f, &jz_f, descZ, work, &lwork, iwork, &liwork, \
                              ifail, iclustr, gap, &info); \
  return info; \
}

CSCALAPACK_PSYEVX_WORK_IMPL(pssyevx, PSSYEVX, float, float)
CSCALAPACK_PSYEVX_WORK_IMPL(pdsyevx, PDSYEVX, double, double)
CSCALAPACK_PSYEVX_WORK_IMPL(pcheevx, PCHEEVX, lapack_complex_float, float)
CSCALAPACK_PSYEVX_WORK_IMPL(pzheevx, PZHEEVX, lapack_complex_double, double)

#undef CSCALAPACK_PSYEVX_WORK_IMPL
