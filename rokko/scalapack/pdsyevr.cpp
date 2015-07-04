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

int ROKKO_pdsyevr(char jobz, char uplo, int n,
		  double* A, int ia, int ja, const int* descA,
		  double* w, double* Z, int iz, int jz, const int* descZ) {

  double* work = new double[1];
  int* iwork = new int[1];
  long lwork = -1;
  long liwork = -1;
  int info;

  // call for querying optimal size of work arrays
  info = ROKKO_pdsyevr_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ,
			    work, lwork, iwork, liwork);

  if (info) {
    std::cerr << "error at querying size. info=" << info  << std::endl;
    exit(1);
  }

  lwork = work[0];
  delete[] work;
  work = new double[lwork];
  liwork = iwork[0];
  delete[] iwork;
  iwork = new int[liwork];
  if (work == 0 || iwork == 0) {
    std::cerr << "failed to allocate work. info=" << info << std::endl;
    return info;
  }

  info = ROKKO_pdsyevd_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ,
			    work, lwork, iwork, liwork);

  if (info) {
    std::cerr << "error at pdsyevd function. info=" << info  << std::endl;
    exit(1);
  }
    
  delete[] work;
  delete[] iwork;

  return info;
}
