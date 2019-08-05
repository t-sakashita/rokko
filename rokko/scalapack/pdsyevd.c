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

int cscalapack_pdsyevd(char jobz, char uplo, int n,
                       double* A, int ia, int ja, const int* descA,
                       double* w, double* Z, int iz, int jz, const int* descZ) {
  // call for querying optimal size of work array
  int lwork = -1;
  int liwork = 1;
  double work_query[1]; 
  int iwork_query[1];
  int info;
  info = cscalapack_pdsyevd_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ,
                                 work_query, lwork, iwork_query, liwork);
  if (info) return info;

  // allocate work arrays
  lwork = (int)work_query[0];
  double* work = (double*)malloc(sizeof(double) * lwork);
  liwork = iwork_query[0];
  int* iwork = (int*)malloc(sizeof(int) * liwork);
  if (work == NULL || iwork == NULL) return 1;

  // call for computation
  info = cscalapack_pdsyevd_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ,
                                 work, lwork, iwork, liwork);
  free(work);
  free(iwork);
  return info;
}
