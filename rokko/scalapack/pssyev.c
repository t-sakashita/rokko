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

int cscalapack_pssyev(char jobz, char uplo, int n,
                      float* A, int ia, int ja, const int* descA,
                      float* w, float* Z, int iz, int jz, const int* descZ) {
  // call for querying optimal size of work array
  int lwork = -1;
  float work_query[1];
  int info;
  info = cscalapack_pssyev_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ,
                                work_query, lwork);
  if (info) return info;

  // allocate work array
  lwork = (int)work_query[0];
  float* work = (float*)malloc(sizeof(float) * lwork);
  if (work == NULL) return 1;

  // call for computation
  info = cscalapack_pssyev_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ,
                                work, lwork);
  free(work);
  return info;
}
