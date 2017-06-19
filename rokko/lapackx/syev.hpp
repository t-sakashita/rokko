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

#ifndef ROKKO_LAPACKX_SYEV_HPP
#define ROKKO_LAPACKX_SYEV_HPP

#include <complex>
#include <lapacke.h>
#include <rokko/traits/norm_t.hpp>
#include <rokko/traits/value_t.hpp>

namespace rokko {
namespace lapackx {

namespace {

template<typename T> struct syev_dispatch;
  
template<>
struct syev_dispatch<float> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int syev(int matrix_layout, char jobz, char uplo, MATRIX& a, VECTOR& w) {
    return LAPACKE_ssyev(matrix_layout, jobz, uplo, rows(a), storage(a), lda(a), &w[0]);
  }
};
  
template<>
struct syev_dispatch<double> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int syev(int matrix_layout, char jobz, char uplo, MATRIX& a, VECTOR& w) {
    return LAPACKE_dsyev(matrix_layout, jobz, uplo, rows(a), storage(a), lda(a), &w[0]);
  }
};
  
template<>
struct syev_dispatch<std::complex<float> > {
  template<typename MATRIX, typename VECTOR>
  static lapack_int syev(int matrix_layout, char jobz, char uplo, MATRIX& a, VECTOR& w) {
    return LAPACKE_cheev(matrix_layout, jobz, uplo, rows(a), storage(a), lda(a), &w[0]);
  }
};
  
template<>
struct syev_dispatch<std::complex<double> > {
  template<typename MATRIX, typename VECTOR>
  static lapack_int syev(int matrix_layout, char jobz, char uplo, MATRIX& a, VECTOR& w) {
    return LAPACKE_zheev(matrix_layout, jobz, uplo, rows(a), storage(a), lda(a), &w[0]);
  }
};

}
  
template<typename MATRIX, typename VECTOR>
lapack_int syev(char jobz, char uplo, MATRIX& a, VECTOR& w) {
  BOOST_STATIC_ASSERT(boost::is_same<typename norm_t<MATRIX>::type, typename value_t<VECTOR>::type>::value);
  return syev_dispatch<typename value_t<MATRIX>::type>::
    syev((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), jobz, uplo, a, w);
}
  
template<typename MATRIX, typename VECTOR>
lapack_int heev(char jobz, char uplo, MATRIX& a, VECTOR& w) {
  BOOST_STATIC_ASSERT(boost::is_same<typename norm_t<MATRIX>::type, typename value_t<VECTOR>::type>::value);
  return syev_dispatch<typename value_t<MATRIX>::type>
    ::syev((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), jobz, uplo, a, w);
}

} // end namespace lapackx
} // end namespace rokko

#endif // ROKKO_LAPACKX_SYEV_HPP
