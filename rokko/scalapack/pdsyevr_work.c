/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/scalapack/scalapack.h>
#include <rokko/scalapack/scalapack_wrap.h>

int ROKKO_pdsyevr_work(char jobz, char range, char uplo, int n,
		       double* A, int ia, int ja, const int* descA,
		       double vl, double vu, int il, int iu,
		       int* m, int* nz,
		       double* w, double* Z, int iz, int jz, const int* descZ,
		       double* work, int lwork, int* iwork, int liwork) {
  int info;
  SCALAPACK_pdsyevr(&jobz, &range, &uplo, &n, A, &ia, &ja, descA, &vl, &vu, &il, &iu, m, nz, w,
                    Z, &iz, &jz, descZ, work, &lwork, iwork, &liwork, &info);
  return info;
}
