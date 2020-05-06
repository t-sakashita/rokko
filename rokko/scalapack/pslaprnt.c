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

void cscalapack_pslaprnt(int m, int n, const float* A, int ia, int ja, const int* descA,
                         int irprnt, int icprnt, const char* cmatnm, int nout, float* work) {
  int ia_f = ia + 1;
  int ja_f = ja + 1;
  SCALAPACK_pslaprnt(&m, &n, A, &ia_f, &ja_f, descA, &irprnt, &icprnt, cmatnm, &nout, work);
}
