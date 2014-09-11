/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/pblas/pblas.h>

void ROKKO_pddot(int N, double* DOT,
                 const double* X, int IX, int JX, int* DESCX, int INCX,
                 const double* Y, int IY, int JY, int* DESCY, int INCY) {
  PBLAS_pddot(&N, DOT, X, &IX, &JX, DESCX, &INCX, Y, &IY, &JY, DESCY, &INCY);
}
