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

#ifndef ROKKO_LAPACKX_HPP
#define ROKKO_LAPACKX_HPP

#include <lapacke.h>

namespace rokko {
namespace lapackx {

template<typename MATRIX, typename VECTOR>
lapack_int getrf(MATRIX& a, VECTOR& ipiv);

template<typename MATRIX0, typename MATRIX1, typename VECTOR>
lapack_int getrs(char trans, lapack_int nrhs, MATRIX0 const& a,
                 VECTOR const& ipiv, MATRIX1& b);

template<typename MATRIX, typename VECTOR>
lapack_int gesvd(char jobu, char jobvt, MATRIX& a, VECTOR& s,
                 MATRIX& u, MATRIX& vt);
  
template<typename MATRIX, typename VECTOR>
lapack_int gesvd_work(char jobu, char jobvt, MATRIX& a, VECTOR& s,
                      MATRIX& u, MATRIX& vt, VECTOR& work);

template<typename MATRIX, typename VECTOR>
lapack_int syev(char jobz, char uplo, MATRIX& a, VECTOR& w);
  
template<typename MATRIX, typename VECTOR>
lapack_int heev(char jobz, char uplo, MATRIX& a, VECTOR& w);

} // end namespace lapackx
} // end namespace rokko

#include "lapackx/gesvd.hpp"
#include "lapackx/getrf.hpp"
#include "lapackx/getrs.hpp"
#include "lapackx/syev.hpp"

#endif // ROKKO_LAPACKX_HPP
