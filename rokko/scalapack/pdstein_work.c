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

int cscalapack_pdstein_work(int n, const double* d, const double* e, int m,
                            double* w, const int* iblock, const int* isplit, double orfac,
                            double* Z, const int* iZ, const int* jZ, const int* descZ,
                            double* work, int lwork, int* iwork, int liwork,
                            int* ifail, int* iclustr, double* gap) {
  int info;
  SCALAPACK_pdstein(&n, d, e, &m,
                    w, iblock, isplit, &orfac,
                    Z, iZ, jZ, descZ,
                    work, &lwork, iwork, &liwork,
                    ifail, iclustr, gap, &info);
  return info;
}
