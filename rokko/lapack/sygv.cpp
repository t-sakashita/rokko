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

#define LAPACKE_SYGV_IMPL(NAME, TYPE) \
lapack_int LAPACKE_sygv(int matrix_order, lapack_int itype, char jobz, char uplo, lapack_int n, TYPE * a, lapack_int lda, TYPE * b, lapack_int ldb, TYPE * w ) { \
  return LAPACKE_ ## NAME (matrix_order, itype, jobz, uplo, n, a, lda, b, ldb, w); }

LAPACKE_SYGV_IMPL(ssygv, float);
LAPACKE_SYGV_IMPL(dsygv, double);

#undef LAPACKE_SYGV_IMPL
