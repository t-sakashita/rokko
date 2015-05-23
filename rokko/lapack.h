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

#ifndef ROKKO_LAPACK_H
#define ROKKO_LAPACK_H

#include <rokko/mangling.h>
#include <lapacke.h>

#ifdef __cplusplus

#include <boost/numeric/ublas/traits.hpp>

#define LAPACKE_GEQRF_DECL(NAME, TYPE) \
lapack_int LAPACKE_geqrf(int matrix_order, lapack_int m, lapack_int n, TYPE * a, lapack_int lda, TYPE * tau);
LAPACKE_GEQRF_DECL(sgeqrf, float);
LAPACKE_GEQRF_DECL(dgeqrf, double);
LAPACKE_GEQRF_DECL(cgeqrf, lapack_complex_float);
LAPACKE_GEQRF_DECL(zgeqrf, lapack_complex_double);
#undef LAPACKE_GEQRF_DECL

#define LAPACKE_GESVD_DECL(NAME, TYPE) \
lapack_int LAPACKE_gesvd(int matrix_order, char jobu, char jobvt, lapack_int m, lapack_int n, TYPE * a, lapack_int lda, boost::numeric::ublas::type_traits<TYPE>::real_type * s, TYPE * u, lapack_int ldu, TYPE * vt, lapack_int ldvt, boost::numeric::ublas::type_traits<TYPE>::real_type * superb);
LAPACKE_GESVD_DECL(sgesvd, float);
LAPACKE_GESVD_DECL(dgesvd, double);
LAPACKE_GESVD_DECL(cgesvd, lapack_complex_float);
LAPACKE_GESVD_DECL(zgesvd, lapack_complex_double);
#undef LAPACKE_GESVD_DECL

#define LAPACKE_HEEV_DECL(NAME, TYPE) \
lapack_int LAPACKE_heev(int matrix_order, char jobz, char uplo, lapack_int n, TYPE * a, lapack_int lda, boost::numeric::ublas::type_traits<TYPE>::real_type * w);
LAPACKE_HEEV_DECL(cheev, lapack_complex_float);
LAPACKE_HEEV_DECL(zheev, lapack_complex_double);
#undef LAPACKE_HEEV_DECL

#define LAPACKE_ORGQR_DECL(NAME, TYPE) \
lapack_int LAPACKE_orgqr(int matrix_order, lapack_int m, lapack_int n, lapack_int k, TYPE * a, lapack_int lda, const TYPE * tau); \
lapack_int LAPACKE_ungqr(int matrix_order, lapack_int m, lapack_int n, lapack_int k, TYPE * a, lapack_int lda, const TYPE * tau);
LAPACKE_ORGQR_DECL(sorgqr, float);
LAPACKE_ORGQR_DECL(dorgqr, double);
#undef LAPACKE_ORGQR_DECL

#define LAPACKE_SYEV_DECL(NAME, TYPE) \
lapack_int LAPACKE_syev(int matrix_order, char jobz, char uplo, lapack_int n, TYPE * a, lapack_int lda, TYPE * w ); \
lapack_int LAPACKE_heev(int matrix_order, char jobz, char uplo, lapack_int n, TYPE * a, lapack_int lda, TYPE * w );
LAPACKE_SYEV_DECL(ssyev, float);
LAPACKE_SYEV_DECL(dsyev, double);
#undef LAPACKE_SYEV_DECL

#define LAPACKE_UNGQR_DECL(NAME, TYPE) \
lapack_int LAPACKE_ungqr(int matrix_order, lapack_int m, lapack_int n, lapack_int k, TYPE * a, lapack_int lda, const TYPE * tau);
LAPACKE_UNGQR_DECL(cungqr, lapack_complex_float);
LAPACKE_UNGQR_DECL(zungqr, lapack_complex_double);
#undef LAPACKE_UNGQR_DECL

#endif // __cplusplus

#endif // ROKKO_LAPACK_H
