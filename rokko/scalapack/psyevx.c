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
#include <stdlib.h>
#include <rokko/scalapack/scalapack_interface.h>

#define CSCALAPACK_PSYEVX_IMPL(NAMES, NAMEL, TYPE, TYPE_REAL) \
int cscalapack_ ## NAMES (char jobz, char range, char uplo, int n, \
                          TYPE* A, int ia, int ja, const int* descA, \
                          TYPE_REAL vl, TYPE_REAL vu, int il, int iu, \
                          TYPE_REAL abstol, int* m, int* nZ, TYPE_REAL* w, TYPE_REAL orfac, \
                          TYPE* Z, int iz, int jz, const int* descZ, \
                          int* ifail, int* iclustr, TYPE_REAL* gap) { \
  /* call for querying optimal size of work array */ \
  int lwork = -1; \
  int liwork = -1; \
  TYPE work_query[1]; \
  int iwork_query[1]; \
  int info; \
  info = cscalapack_## NAMES ##_work(jobz, range, uplo, n, A, ia, ja, descA, \
                                     vl, vu, il, iu, \
                                     abstol, m, nZ, w, orfac, \
                                     Z, iz, jz, descZ, \
                                     work_query, lwork, iwork_query, liwork, \
                                     ifail, iclustr, gap); \
  if (info) return info; \
  \
  /* allocate work arrays */ \
  lwork = (int)work_query[0]; \
  TYPE* work = (TYPE*)malloc( sizeof(TYPE) * lwork ); \
  liwork = iwork_query[0]; \
  int* iwork = (int*)malloc( sizeof(int) * liwork ); \
  if (work == NULL || iwork == NULL) return 1; \
  \
  /* call for computation */ \
  info = cscalapack_## NAMES ##_work(jobz, range, uplo, n, A, ia, ja, descA, \
                                    vl, vu, il, iu, \
                                    abstol, m, nZ, w, orfac, \
                                    Z, iz, jz, descZ, \
                                    work, lwork, iwork, liwork, \
                                    ifail, iclustr, gap); \
  free(work); \
  free(iwork); \
  return info; \
}

CSCALAPACK_PSYEVX_IMPL(pssyevx, PSSYEVX, float, float)
CSCALAPACK_PSYEVX_IMPL(pdsyevx, PDSYEVX, double, double)
CSCALAPACK_PSYEVX_IMPL(pcheevx, PCHEEVX, lapack_complex_float, float)
CSCALAPACK_PSYEVX_IMPL(pzheevx, PZHEEVX, lapack_complex_double, double)

#undef CSCALAPACK_PSYEVX_IMPL
