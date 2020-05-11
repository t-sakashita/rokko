/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/cscalapack.h>
#include <stdlib.h>
#include <rokko/scalapack/scalapack_interface.h>

#define CSCALAPACK_PSYEVD_IMPL(NAMES, NAMEL, TYPE) \
int cscalapack_ ## NAMES (char jobz, char uplo, int n, \
                          TYPE* A, int ia, int ja, const int* descA,  \
                          TYPE* w, TYPE* Z, int iz, int jz, const int* descZ) { \
  /* call for querying optimal size of work array */                    \
  int lwork = -1; \
  int liwork = 1; \
  TYPE work_query[1]; \
  int iwork_query[1]; \
  int info; \
  info = cscalapack_ ## NAMES ## _work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ, \
                                 work_query, lwork, iwork_query, liwork); \
  if (info) return info; \
  \
  /* allocate work arrays */  \
  lwork = (int)work_query[0]; \
  TYPE* work = (TYPE*)malloc(sizeof(TYPE) * lwork); \
  liwork = iwork_query[0]; \
  int* iwork = (int*)malloc(sizeof(int) * liwork); \
  if (work == NULL || iwork == NULL) return 1; \
  \
  /* call for computation */                                             \
  info = cscalapack_ ## NAMES ## _work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ, \
                                 work, lwork, iwork, liwork); \
  free(work); \
  free(iwork); \
  return info; \
}

CSCALAPACK_PSYEVD_IMPL(pssyevd, PSSYEVD, float)
CSCALAPACK_PSYEVD_IMPL(pdsyevd, PDSYEVD, double)
CSCALAPACK_PSYEVD_IMPL(pcheevd, PCHEEVD, lapack_complex_float)
CSCALAPACK_PSYEVD_IMPL(pzheevd, PZHEEVD, lapack_complex_double)

#undef CSCALAPACK_PSYEVD_IMPL
