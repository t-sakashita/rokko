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

double cscalapack_pdlange_work(char norm, int m, int n, const double* A, const int* descA,
                               double* work) {
  int ia = 1;
  int ja = 1;
  return SCALAPACK_pdlange(&norm, &m, &n, A, &ia, &ja, descA, work);
}
