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

#ifndef ROKKO_LAPACK_GEQRF_HPP
#define ROKKO_LAPACK_GEQRF_HPP

#include <complex>
#include <stdexcept>
#include <lapacke.h>
#include "complex_cast.hpp"

namespace rokko {
namespace lapack {

namespace {

template<typename T> struct geqrf_dispatch;
  
template<>
struct geqrf_dispatch<float> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int geqrf(int matrix_layout, lapack_int m, lapack_int n,
                          MATRIX& a, VECTOR& tau) {
    return LAPACKE_sgeqrf(matrix_layout, m, n, storage(a), lda(a), storage(tau));
  }
  template<typename MATRIX, typename VECTOR>
  static lapack_int geqrf(int matrix_layout, lapack_int m, lapack_int n,
                          MATRIX& a, VECTOR& tau, VECTOR& work) {
    return LAPACKE_sgeqrf(matrix_layout, m, n, storage(a), lda(a), storage(tau),
                          storage(work), size(work));
  }
};

template<>
struct geqrf_dispatch<double> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int geqrf(int matrix_layout, lapack_int m, lapack_int n,
                          MATRIX& a, VECTOR& tau) {
    return LAPACKE_dgeqrf(matrix_layout, m, n, storage(a), lda(a), storage(tau));
  }
  template<typename MATRIX, typename VECTOR>
  static lapack_int geqrf(int matrix_layout, lapack_int m, lapack_int n,
                          MATRIX& a, VECTOR& tau, VECTOR& work) {
    return LAPACKE_dgeqrf(matrix_layout, m, n, storage(a), lda(a), storage(tau),
                          storage(work), size(work));
  }
};

template<>
struct geqrf_dispatch<std::complex<float> > {
  template<typename MATRIX, typename VECTOR>
  static lapack_int geqrf(int matrix_layout, lapack_int m, lapack_int n,
                          MATRIX& a, VECTOR& tau) {
    return LAPACKE_cgeqrf(matrix_layout, m, n, complex_cast(storage(a)), lda(a),
                          complex_cast(storage(tau)));
  }
  template<typename MATRIX, typename VECTOR>
  static lapack_int geqrf(int matrix_layout, lapack_int m, lapack_int n,
                          MATRIX& a, VECTOR& tau, VECTOR& work) {
    return LAPACKE_cgeqrf(matrix_layout, m, n, complex_cast(storage(a)), lda(a),
                          complex_cast(storage(tau)), complex_cast(storage(work)),
                          size(work));
  }
};

template<>
struct geqrf_dispatch<std::complex<double> > {
  template<typename MATRIX, typename VECTOR>
  static lapack_int geqrf(int matrix_layout, lapack_int m, lapack_int n,
                          MATRIX& a, VECTOR& tau) {
    return LAPACKE_zgeqrf(matrix_layout, m, n, complex_cast(storage(a)), lda(a),
                          complex_cast(storage(tau)));
  }
  template<typename MATRIX, typename VECTOR>
  static lapack_int geqrf(int matrix_layout, lapack_int m, lapack_int n,
                          MATRIX& a, VECTOR& tau, VECTOR& work) {
    return LAPACKE_zgeqrf(matrix_layout, m, n, complex_cast(storage(a)), lda(a),
                          complex_cast(storage(tau)), complex_cast(storage(work)),
                          size(work));
  }
};

}
  
template<typename MATRIX, typename VECTOR>
lapack_int geqrf(MATRIX& a, VECTOR& tau) {
  BOOST_STATIC_ASSERT(boost::is_same<typename value_t<MATRIX>::type,
                      typename value_t<VECTOR>::type>::value);
  lapack_int m = rows(a);
  lapack_int n = cols(a);
  lapack_int r = std::min(m, n);
  if (size(tau) != r)
    throw std::invalid_argument("vector tau size mismatch");
  return geqrf_dispatch<typename value_t<MATRIX>::type>
    ::geqrf((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), m, n, a, tau);
}

template<typename MATRIX, typename VECTOR>
lapack_int geqrf(MATRIX& a, VECTOR& tau, VECTOR& work) {
  BOOST_STATIC_ASSERT(boost::is_same<typename value_t<MATRIX>::type,
                      typename value_t<VECTOR>::type>::value);
  lapack_int m = rows(a);
  lapack_int n = cols(a);
  lapack_int r = std::min(m, n);
  if (size(tau) != r)
    throw std::invalid_argument("vector tau size mismatch");
  if (size(work) < std::max(1, n))
    throw std::invalid_argument("vector work size mismatch");
  return geqrf_dispatch<typename value_t<MATRIX>::type>
    ::geqrf((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), m, n, a, tau, work);
}

} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_GEQRF_HPP
