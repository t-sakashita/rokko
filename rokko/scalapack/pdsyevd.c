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

#include <stdio.h>
#include <stdlib.h>
#include <rokko/scalapack/scalapack.h>
#include <rokko/scalapack/scalapack_wrap.h>

int ROKKO_pdsyevd(char jobz, char uplo, int n,
		  double* A, int ia, int ja, const int* descA,
		  double* w, double* Z, int iz, int jz, const int* descZ) {
  // call for querying optimal size of work array
  int lwork = -1;
  int liwork = 1;
  double work_query[1]; 
  int iwork_query[1];
  int info;
  info = ROKKO_pdsyevd_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ,
			    work_query, lwork, iwork_query, liwork);
  if (info) {
    printf("error at querying size at pdsyevd. info=%d\n", info);
    exit(1);
  }

  // allocate work arrays
  lwork = (int)work_query[0];
  double* work = (double*)malloc( sizeof(double) * lwork );
  liwork = iwork_query[0];
  int* iwork = (int*)malloc( sizeof(int) * liwork );
  if (work == NULL || iwork == NULL) {
    printf("failed to allocate work. info=%d\n", info);
    return info;
  }

  // call for computation
  info = ROKKO_pdsyevd_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ,
			    work, lwork, iwork, liwork);
  if (info) {
    printf("error at pdsyevd function. info=%d\n", info);
    exit(1);
  }
    
  free(work);
  free(iwork);
  return info;
}
