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

#define LAPACKE_SYEV_IMPL(NAME, TYPE) \
lapack_int LAPACKE_heev(int matrix_order, char jobz, char uplo, lapack_int n, TYPE * a, lapack_int lda, TYPE * w ) { \
  return LAPACKE_ ## NAME (matrix_order, jobz, uplo, n, a, lda, w ); }

LAPACKE_SYEV_IMPL(ssyev, float);
LAPACKE_SYEV_IMPL(dsyev, double);

#undef LAPACKE_SYEV_IMPL

#define LAPACKE_HEEV_IMPL(NAME, TYPE) \
lapack_int LAPACKE_heev(int matrix_order, char jobz, char uplo, lapack_int n, TYPE * a, lapack_int lda, boost::numeric::ublas::type_traits<TYPE>::real_type * w) { \
  return LAPACKE_ ## NAME (matrix_order, jobz, uplo, n, a, lda, w); }

LAPACKE_HEEV_IMPL(cheev, lapack_complex_float);
LAPACKE_HEEV_IMPL(zheev, lapack_complex_double);

#undef LAPACKE_HEEV_IMPL
