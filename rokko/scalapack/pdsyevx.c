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

int ROKKO_pdsyevx(char jobz, char range, char uplo, int n,
		  double* A, int ia, int ja, const int* descA,
		  double vl, double vu, int il, int iu,
		  double abstol, int* m, int* nZ, double* w, double orfac,
		  double* Z, int iz, int jz, const int* descZ,
		  int* ifail, int* iclustr, double* gap) {
  // call for querying optimal size of work array
  int lwork = -1;
  int liwork = -1;
  double work_query[1]; 
  int iwork_query[1];
  int info;
  info = ROKKO_pdsyevx_work(jobz, range, uplo, n, A, ia, ja, descA,
			    vl, vu, il, iu,
			    abstol, m, nZ, w, orfac,
			    Z, iz, jz, descZ,
			    work_query, lwork, iwork_query, liwork,
			    ifail, iclustr, gap);
  if (info) {
    printf("error at querying size at pdsyevdr. info=%d\n", info);
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
  info = ROKKO_pdsyevx_work(jobz, range, uplo, n, A, ia, ja, descA,
			    vl, vu, il, iu,
			    abstol, m, nZ, w, orfac,
			    Z, iz, jz, descZ,
			    work, lwork, iwork, liwork,
			    ifail, iclustr, gap);
  printf("reeeeeached info=%d\n", info);

  if (info) {
    printf("error at pdsyevx function. info=%d\n", info);
    exit(1);
  }
    
  free(work);
  free(iwork);
  printf("qqqqqqqqqqq\n");
  return info;
}
