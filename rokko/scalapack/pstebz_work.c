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

#define CSCALAPACK_PSTEBZ_WORK_IMPL(NAMES, NAMEL, TYPE) \
int cscalapack_## NAMES ##_work(int ictxt, char range, char order, int n, \
                                TYPE vl, TYPE vu, int il, int iu, \
                                TYPE abstol, const TYPE* d, const TYPE* e, int* m, int* nsplit, \
                                TYPE* w, int* iblock, int* isplit, \
                                TYPE* work, int lwork, int* iwork, int liwork) { \
  int info; \
  ROKKO_GLOBAL(NAMES, NAMEL) (&ictxt, &range, &order, &n, \
                              &vl, &vu, &il, &iu, \
                              &abstol, d, e, m, nsplit, \
                              w, iblock, isplit, \
                              work, &lwork, iwork, &liwork, &info); \
  return info; \
}

CSCALAPACK_PSTEBZ_WORK_IMPL(psstebz, PSSTEBZ, float)
CSCALAPACK_PSTEBZ_WORK_IMPL(pdstebz, PDSTEBZ, double)

#undef CSCALAPACK_PSTEBZ_WORK_IMPL
