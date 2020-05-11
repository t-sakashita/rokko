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

#define CSCALAPACK_PSYEVR_WORK_IMPL(NAMES, NAMEL, TYPE, TYPE_REAL) \
int cscalapack_## NAMES ##_work(char jobz, char range, char uplo, int n, \
                                TYPE* A, int ia, int ja, const int* descA, \
                                TYPE_REAL vl, TYPE_REAL vu, int il, int iu, \
                                int* m, int* nz, \
                                TYPE_REAL* w, TYPE* Z, int iz, int jz, const int* descZ, \
                                TYPE* work, int lwork, int* iwork, int liwork) { \
  int ia_f = ia + 1; \
  int ja_f = ja + 1; \
  int il_f = il + 1; \
  int iu_f = iu + 1; \
  int iz_f = iz + 1; \
  int jz_f = jz + 1; \
  int info; \
  ROKKO_GLOBAL(NAMES, NAMEL) (&jobz, &range, &uplo, &n, A, &ia_f, &ja_f, descA, &vl, &vu, &il_f, &iu_f, \
                              m, nz, w, Z, &iz_f, &jz_f, descZ, work, &lwork, iwork, &liwork, &info); \
  return info; \
}

CSCALAPACK_PSYEVR_WORK_IMPL(pssyevr, PSSYEVR, float, float)
CSCALAPACK_PSYEVR_WORK_IMPL(pdsyevr, PDSYEVR, double, double)
CSCALAPACK_PSYEVR_WORK_IMPL(pcheevr, PCHEEVR, lapack_complex_float, float)
CSCALAPACK_PSYEVR_WORK_IMPL(pzheevr, PZHEEVR, lapack_complex_double, double)

#undef CSCALAPACK_PSYEVR_WORK_IMPL
