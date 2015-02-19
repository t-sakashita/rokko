/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2015 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/pblas/pblas.h>
#include <rokko/pblas/pblas_wrap.h>

void ROKKO_pdcopy(int N, const double* X, int IX, int JX, int* DESCX, int INCX,
                  double* Y, int IY, int JY, int *DESCY, int INCY) {
  PBLAS_pdcopy(&N, X, &IX, &JX, DESCX, &INCX, Y, &IY, &JY, DESCY, &INCY);
}
