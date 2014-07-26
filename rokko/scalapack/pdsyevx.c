/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/scalapack/scalapack.h>

void ROKKO_pdsyevx(char jobz, char range, char uplo, int n,
                   double* A, int iA, int jA, const int* descA,
                   double vl, double vu, int il, int iu,
                   double abstol, int* m, int* nZ, double* w, double orfac,
                   double* Z, int iZ, int jZ, const int* descZ,
                   double* work, int lwork, int* iwork, int liwork,
                   int* ifail, int* iclustr, double* gap, int* info) {
  pdsyevx_(&jobz, &range, &uplo, &n, A, &iA, &jA, descA, &vl, &vu, &il, &iu,
           &abstol, m, nZ, w, &orfac, Z, &iZ, &jZ, descZ, work, &lwork, iwork, &liwork,
           ifail, iclustr, gap, info);
}
