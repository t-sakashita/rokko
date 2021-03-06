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

void cscalapack_pselset(float* A, int ia, int ja, const int* descA, float alpha) {
  int ia_f = ia + 1;
  int ja_f = ja + 1;
  SCALAPACK_pselset(A, &ia_f, &ja_f, descA, &alpha);
}
