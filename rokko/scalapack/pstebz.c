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

#define CSCALAPACK_PSTEBZ_IMPL(NAMES, NAMEL, TYPE) \
int cscalapack_## NAMES (int ictxt, char range, char order, int n, \
                         TYPE vl, TYPE vu, int il, int iu, \
                         TYPE abstol, const TYPE* d, const TYPE* e, int* m, int* nsplit, \
                         TYPE* w, int* iblock, int* isplit) { \
  /* call for querying optimal size of work array */ \
  int lwork = -1; \
  int liwork = -1; \
  TYPE work_query[1]; \
  int iwork_query[1]; \
  int info; \
  info = cscalapack_## NAMES ##_work(ictxt, range, order, n, \
                                     vl, vu, il, iu, \
                                     abstol, d, e, m, nsplit, \
                                     w, iblock, isplit, \
                                     work_query, lwork, iwork_query, liwork); \
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
  info = cscalapack_## NAMES ##_work(ictxt, range, order, n, \
                                     vl, vu, il, iu, \
                                     abstol, d, e, m, nsplit, \
                                     w, iblock, isplit, \
                                     work, lwork, iwork, liwork); \
  free(work); \
  free(iwork); \
  return info; \
}

CSCALAPACK_PSTEBZ_IMPL(psstebz, PSSTEBZ, float)
CSCALAPACK_PSTEBZ_IMPL(pdstebz, PDSTEBZ, double)

#undef CSCALAPACK_PSTEBZ_IMPL
