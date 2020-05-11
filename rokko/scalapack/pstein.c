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

#define CSCALAPACK_PSTEIN_IMPL(NAMES, NAMEL, TYPE, TYPE_REAL) \
int cscalapack_## NAMES (int n, const TYPE_REAL* d, const TYPE_REAL* e, int m, \
                         TYPE_REAL* w, const int* iblock, const int* isplit, TYPE_REAL orfac, \
                         TYPE* Z, const int* iZ, const int* jZ, const int* descZ, \
                         int* ifail, int* iclustr, TYPE_REAL* gap) { \
   /* call for querying optimal size of work array */ \
  int lwork = -1; \
  int liwork = -1; \
  TYPE work_query[1]; \
  int iwork_query[1]; \
  int info; \
  info = cscalapack_## NAMES ##_work(n, d, e, m, \
                                     w, iblock, isplit, orfac, \
                                     Z, iZ, jZ, descZ, \
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
  info = cscalapack_## NAMES ##_work(n, d, e, m, \
                                     w, iblock, isplit, orfac, \
                                     Z, iZ, jZ, descZ, \
                                     work, lwork, iwork, liwork, \
                                     ifail, iclustr, gap); \
  free(work); \
  free(iwork); \
  return info; \
}

CSCALAPACK_PSTEIN_IMPL(psstein, PSSTEIN, float, float)
CSCALAPACK_PSTEIN_IMPL(pdstein, PDSTEIN, double, double)
CSCALAPACK_PSTEIN_IMPL(pcstein, PCSTEIN, lapack_complex_float, float)
CSCALAPACK_PSTEIN_IMPL(pzstein, PZSTEIN, lapack_complex_double, double)

#undef CSCALAPACK_PSTEIN_IMPL
