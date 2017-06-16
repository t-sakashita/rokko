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
#include <boost/numeric/ublas/traits.hpp>
#include "matrix_traits.hpp"

namespace rokko {
namespace lapackx {

template<typename MATRIX, typename VECTOR>
lapack_int getrf(MATRIX& a, VECTOR& ipiv);

template<typename MATRIX0, typename VECTOR, typename MATRIX1>
lapack_int getrs(char trans, lapack_int nrhs, MATRIX0 const& a,
                 VECTOR const& ipiv, MATRIX1& b);

// template<typename MATRIX, typename VECTOR>
// lapack_int gesvd(char jobu, char jobvt, MATRIX& a, VECTOR& s,
//                  MATRIX& u, MATRIX& vt);
  
template<typename MATRIX, typename VECTOR>
lapack_int gesvd_work(char jobu, char jobvt, MATRIX& a, VECTOR& s,
                      MATRIX& u, MATRIX& vt, VECTOR& work);

template<typename MATRIX>
typename boost::numeric::ublas::type_traits<
  typename matrix_traits<MATRIX>::value_type>::real_type
lange(char norm, MATRIX const& a);
  
template<typename MATRIX, typename VECTOR>
typename boost::numeric::ublas::type_traits<
  typename matrix_traits<MATRIX>::value_type>::real_type
lange_work(char norm, MATRIX const& a, VECTOR& work);
  
template<typename MATRIX, typename VECTOR>
lapack_int syev(char jobz, char uplo, MATRIX& a, VECTOR& w);
  
template<typename MATRIX, typename VECTOR>
lapack_int heev(char jobz, char uplo, MATRIX& a, VECTOR& w);

} // end namespace lapackx
} // end namespace rokko

#include "lapackx/gesvd.hpp"
#include "lapackx/syev.hpp"

#endif // ROKKO_LAPACKX_HPP
