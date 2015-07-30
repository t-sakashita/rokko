/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2015 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/lapack.h>

#define LAPACKE_GEQRF_IMPL(NAME, TYPE) \
lapack_int LAPACKE_geqrf(int matrix_order, lapack_int m, lapack_int n, TYPE * a, lapack_int lda, TYPE * tau) { \
  return LAPACKE_ ## NAME (matrix_order, m, n, a, lda, tau); }

LAPACKE_GEQRF_IMPL(sgeqrf, float);
LAPACKE_GEQRF_IMPL(dgeqrf, double);
LAPACKE_GEQRF_IMPL(cgeqrf, lapack_complex_float);
LAPACKE_GEQRF_IMPL(zgeqrf, lapack_complex_double);

#undef LAPACKE_GEQRF_IMPL
