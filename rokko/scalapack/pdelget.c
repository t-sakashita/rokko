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

void ROKKO_pdelget(char scope, char top, double* alpha,
                   const double* A, int ia, int ja, const int* descA) {
  SCALAPACK_pdelget(&scope, &top, alpha, A, &ia, &ja, descA);
}
