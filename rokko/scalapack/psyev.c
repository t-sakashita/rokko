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

#define CSCALAPACK_PSYEV_IMPL(NAMES, NAMEL, TYPE) \
int cscalapack_ ## NAMES (char jobz, char uplo, int n, \
                          TYPE* A, int ia, int ja, const int* descA, \
                          TYPE* w, TYPE* Z, int iz, int jz, const int* descZ) { \
  /* call for querying optimal size of work array */ \
  int lwork = -1; \
  TYPE work_query[1]; \
  int info; \
  info = cscalapack_## NAMES ##_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ, \
                                     work_query, lwork); \
  if (info) return info; \
  \
  /* allocate work array */ \
  lwork = (int)work_query[0]; \
  TYPE* work = (TYPE*)malloc(sizeof(TYPE) * lwork); \
  if (work == NULL) return 1; \
  \
  /* call for computation */ \
  info = cscalapack_ ## NAMES ## _work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ, \
                                       work, lwork); \
  free(work); \
  return info; \
}

CSCALAPACK_PSYEV_IMPL(pssyev, PSSYEV, float)
CSCALAPACK_PSYEV_IMPL(pdsyev, PDSYEV, double)
CSCALAPACK_PSYEV_IMPL(pcheev, PCHEEV, lapack_complex_float)
CSCALAPACK_PSYEV_IMPL(pzheev, PZHEEV, lapack_complex_double)

#undef CSCALAPACK_PSYEV_IMPL
