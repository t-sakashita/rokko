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

#define LAPACKE_GESVD_IMPL(NAME, TYPE) \
lapack_int LAPACKE_gesvd(int matrix_order, char jobu, char jobvt, lapack_int m, lapack_int n, TYPE * a, lapack_int lda, boost::numeric::ublas::type_traits<TYPE>::real_type * s, TYPE * u, lapack_int ldu, TYPE * vt, lapack_int ldvt, boost::numeric::ublas::type_traits<TYPE>::real_type * superb) { \
  return LAPACKE_ ## NAME (matrix_order, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb); }

LAPACKE_GESVD_IMPL(sgesvd, float);
LAPACKE_GESVD_IMPL(dgesvd, double);
LAPACKE_GESVD_IMPL(cgesvd, lapack_complex_float);
LAPACKE_GESVD_IMPL(zgesvd, lapack_complex_double);

#undef LAPACKE_GESVD_IMPL
