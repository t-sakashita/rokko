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
#include <rokko/scalapack/scalapack_interface.h>

int cscalapack_pdsyevx_work(char jobz, char range, char uplo, int n,
                            double* A, int ia, int ja, const int* descA,
                            double vl, double vu, int il, int iu,
                            double abstol, int* m, int* nZ, double* w, double orfac,
                            double* Z, int iz, int jz, const int* descZ,
                            double* work, int lwork, int* iwork, int liwork,
                            int* ifail, int* iclustr, double* gap) {
  int ia_f = ia + 1;
  int ja_f = ja + 1;
  int il_f = il + 1;
  int iu_f = iu + 1;
  int iz_f = iz + 1;
  int jz_f = jz + 1;
  int info;
  SCALAPACK_pdsyevx(&jobz, &range, &uplo, &n, A, &ia_f, &ja_f, descA, &vl, &vu, &il_f, &iu_f,
                    &abstol, m, nZ, w, &orfac, Z, &iz_f, &jz_f, descZ, work, &lwork, iwork, &liwork,
                    ifail, iclustr, gap, &info);
  return info;
}
