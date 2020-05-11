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

#define CSCALAPACK_PSTEIN_WORK_IMPL(NAMES, NAMEL, TYPE, TYPE_REAL) \
int cscalapack_## NAMES ##_work(int n, const TYPE_REAL* d, const TYPE_REAL* e, int m, \
                                TYPE_REAL* w, const int* iblock, const int* isplit, TYPE_REAL orfac, \
                                TYPE* Z, const int* iZ, const int* jZ, const int* descZ, \
                                TYPE* work, int lwork, int* iwork, int liwork, \
                                int* ifail, int* iclustr, TYPE_REAL* gap) { \
  int info; \
  ROKKO_GLOBAL(NAMES, NAMEL) (&n, d, e, &m, \
                              w, iblock, isplit, &orfac, \
                              Z, iZ, jZ, descZ, \
                              work, &lwork, iwork, &liwork, \
                              ifail, iclustr, gap, &info); \
  return info; \
}

CSCALAPACK_PSTEIN_WORK_IMPL(psstein, PSSTEIN, float, float)
CSCALAPACK_PSTEIN_WORK_IMPL(pdstein, PDSTEIN, double, double)
CSCALAPACK_PSTEIN_WORK_IMPL(pcstein, PCSTEIN, lapack_complex_float, float)
CSCALAPACK_PSTEIN_WORK_IMPL(pzstein, PZSTEIN, lapack_complex_double, double)

#undef CSCALAPACK_PSTEIN_WORK_IMPL
