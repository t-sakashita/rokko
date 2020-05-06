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
#include <rokko/scalapack/scalapack_interface.h>

int cscalapack_pssyevd_work(char jobz, char uplo, int n,
                            float* A, int ia, int ja, const int* descA,
                            float* w, float* Z, int iz, int jz, const int* descZ,
                            float* work, int lwork, int* iwork, int liwork) {
  int ia_f = ia + 1;
  int ja_f = ja + 1;
  int iz_f = iz + 1;
  int jz_f = jz + 1;
  int info;
  SCALAPACK_pssyevd(&jobz, &uplo, &n, A, &ia_f, &ja_f, descA, w, Z, &iz_f, &jz_f, descZ,
                    work, &lwork, iwork, &liwork, &info);
  return info;
}
