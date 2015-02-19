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

#include <rokko/scalapack/scalapack.h>
#include <rokko/scalapack/scalapack_wrap.h>

int ROKKO_pdsyev(char jobz, char uplo, int n,
                 double* A, int ia, int ja, const int* descA,
                 double* w, double* Z, int iz, int jz, const int* descZ,
                 double* work, int lwork) {
  int info;
  SCALAPACK_pdsyev(&jobz, &uplo, &n, A, &ia, &ja, descA, w, Z, &iz, &jz, descZ,
                   work, &lwork, &info);
  return info;
}
