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

int cscalapack_pdstebz_work(int ictxt, char range, char order, int n,
                            double vl, double vu, int il, int iu,
                            double abstol, const double* d, const double* e, int* m, int* nsplit,
                            double* w, int* iblock, int* isplit,
                            double* work, int lwork, int* iwork, int liwork) {
  int info;
  SCALAPACK_pdstebz(&ictxt, &range, &order, &n,
                    &vl, &vu, &il, &iu,
                    &abstol, d, e, m, nsplit,
                    w, iblock, isplit,
                    work, &lwork, iwork, &liwork, &info);
  return info;
}
