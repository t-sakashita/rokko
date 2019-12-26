/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2017 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_LAPACK_HPP
#define ROKKO_LAPACK_HPP

#include <lapacke.h>
#undef I
#include <rokko/traits/norm_t.hpp>

namespace rokko {
namespace lapack {

// LAPACKE_*gesvd is called if size(work) == min(rows(a),cols(a))-1,
// otherwise LAPACKE_sgesvd_work or LAPACKE_dgesvd_work is called
template<typename MATRIX, typename VECTOR>
lapack_int gesvd(char jobu, char jobvt, MATRIX& a, VECTOR& s, MATRIX& u, MATRIX& vt,
                 VECTOR& work);

// LAPACKE_cgesvd_work or LAPACKE_zgesvd_work is called
template<typename MATRIX, typename VECTOR0, typename VECTOR1>
lapack_int gesvd(char jobu, char jobvt, MATRIX& a, VECTOR0& s, MATRIX& u, MATRIX& vt,
                 VECTOR1& work, VECTOR0& rwork);

template<typename MATRIX, typename VECTOR>
lapack_int getrf(MATRIX& a, VECTOR& ipiv);

template<typename MATRIX, typename VECTOR>
lapack_int dgetri(MATRIX& a, VECTOR const& ipiv);

template<typename MATRIX, typename VECTOR0, typename VECTOR1>
lapack_int dgetri(MATRIX& a, VECTOR0 const& ipiv, VECTOR1& work);

template<typename MATRIX0, typename MATRIX1, typename VECTOR>
lapack_int getrs(char trans, lapack_int nrhs, MATRIX0 const& a, VECTOR const& ipiv,
                 MATRIX1& b);

template<typename MATRIX, typename VECTOR>
lapack_int geqrf(MATRIX& a, VECTOR& tau);

template<typename MATRIX, typename VECTOR>
lapack_int geqrf(MATRIX& a, VECTOR& tau, VECTOR& work);
  
template<typename MATRIX>
rokko::norm_t<MATRIX> lange(char norm, MATRIX const& a);

template<typename MATRIX, typename VECTOR>
rokko::norm_t<MATRIX> lange(char norm, MATRIX const& a, VECTOR& work);
  
template<typename MATRIX, typename VECTOR>
lapack_int orgqr(lapack_int k, MATRIX& a, VECTOR const& tau);
  
template<typename MATRIX, typename VECTOR>
lapack_int orgqr(lapack_int k, MATRIX& a, VECTOR const& tau, VECTOR& work);
  
template<typename MATRIX, typename VECTOR>
lapack_int ungqr(lapack_int k, MATRIX& a, VECTOR const& tau);
  
template<typename MATRIX, typename VECTOR>
lapack_int ungqr(lapack_int k, MATRIX& a, VECTOR const& tau, VECTOR& work);
  
template<typename MATRIX, typename VECTOR>
lapack_int heev(char jobz, char uplo, MATRIX& a, VECTOR& w);

template<typename MATRIX, typename VECTOR0, typename VECTOR1>
lapack_int heev(char jobz, char uplo, MATRIX& a, VECTOR0& w, VECTOR1& work, VECTOR1& rwork);
  
template<typename MATRIX, typename VECTOR>
lapack_int hegv(lapack_int itype, char jobz, char uplo, MATRIX& a, VECTOR& b, VECTOR& w);
  
template<typename MATRIX, typename VECTOR0, typename VECTOR1>
lapack_int hegv(lapack_int itype, char jobz, char uplo, MATRIX& a, VECTOR0& b, VECTOR0& w,
                VECTOR1& work);

template<typename MATRIX, typename VECTOR>
lapack_int syev(char jobz, char uplo, MATRIX& a, VECTOR& w);

template<typename MATRIX, typename VECTOR0, typename VECTOR1>
lapack_int syev(char jobz, char uplo, MATRIX& a, VECTOR0& w, VECTOR1& work, VECTOR1& rwork);

template<typename MATRIX, typename VECTOR>
lapack_int sygv(lapack_int itype, char jobz, char uplo, MATRIX& a, VECTOR& b,
                VECTOR& w);
  
template<typename MATRIX, typename VECTOR0, typename VECTOR1>
lapack_int sygv(lapack_int itype, char jobz, char uplo, MATRIX& a, VECTOR0& b,
                VECTOR0& w, VECTOR1& work);
  
} // end namespace lapack
} // end namespace rokko

#include "lapack/geqrf.hpp"
#include "lapack/gesvd.hpp"
#include "lapack/getrf.hpp"
#include "lapack/getrs.hpp"
#include "lapack/heev.hpp"
// #include "lapack/hegv.hpp"
#include "lapack/lange.hpp"
#include "lapack/ungqr.hpp"

#endif // ROKKO_LAPACK_HPP
