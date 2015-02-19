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

void ROKKO_pdlaprnt(int m, int n, const double* A, int ia, int ja, const int* descA, int irprnt,
                    int icprnt, const char* cmatnm, int nout, double* work) {
  SCALAPACK_pdlaprnt(&m, &n, A, &ia, &ja, descA, &irprnt, &icprnt, cmatnm, &nout, work);
}
