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

#include <iostream>
#include <rokko/scalapack/scalapack.h>
#include <rokko/scalapack/scalapack_wrap.h>

int ROKKO_pdsyev(char jobz, char uplo, int n,
		 double* A, int ia, int ja, const int* descA,
		 double* w, double* Z, int iz, int jz, const int* descZ) {

  // call for querying optimal size of work array
  double* work = new double[1];
  long lwork = -1;
  int info;
  info = ROKKO_pdsyev_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ,
			   work, lwork);
  if (info) {
    std::cerr << "error at querying size. info=" << info  << std::endl;
    exit(1);
  }

  // allocate work array
  lwork = work[0];
  delete[] work;
  work = new double[lwork];
  if (work == 0) {
    std::cerr << "failed to allocate work. info=" << info << std::endl;
    return info;
  }

  // call for computation
  info = ROKKO_pdsyev_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ,
			   work, lwork);
  if (info) {
    std::cerr << "error at pdsyevd function. info=" << info  << std::endl;
    exit(1);
  }
    
  delete[] work;

  return info;
}
