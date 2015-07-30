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

int ROKKO_pdsyev(char jobz, char uplo, int n,
		 double* A, int ia, int ja, const int* descA,
		 double* w, double* Z, int iz, int jz, const int* descZ) {

  // call for querying optimal size of work array
  int lwork = -1;
  double work_query[1]; 
  int info;
  info = ROKKO_pdsyev_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ,
			   work_query, lwork);
  if (info) {
    printf("error at querying size at pdsyev. info=%d\n", info);
    exit(1);
  }

  // allocate work array
  lwork = (int)work_query[0];
  double* work = (double*)malloc( sizeof(double) * lwork );
  if (work == NULL) {
    printf("failed to allocate work. info=%d\n", info);
    return info;
  }
  
  // call for computation
  info = ROKKO_pdsyev_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ,
			   work, lwork);
  if (info) {
    printf("error at pdsyev function. info=%d\n", info);
    exit(1);
  }
    
  free(work);
  return info;
}
