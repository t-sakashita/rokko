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

#ifndef ROKKO_LAPACKX_GETRF_HPP
#define ROKKO_LAPACKX_GETRF_HPP

#include <rokko/vector_traits.hpp>
#include <rokko/matrix_traits.hpp>
#include <complex>
#include <stdexcept>
#include <lapacke.h>

namespace rokko {
namespace lapackx {

namespace {

template<typename T1, typename T2> struct getrf_dispatch;
  
template<>
struct getrf_dispatch<float, lapack_int> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int getrf(int matrix_layout, lapack_int m, lapack_int n, MATRIX& a, VECTOR& ipiv) {
    return LAPACKE_sgetrf(matrix_layout, m, n, storage(a), lda(a), storage(ipiv));
  }
};

template<>
struct getrf_dispatch<double, lapack_int> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int getrf(int matrix_layout, lapack_int m, lapack_int n, MATRIX& a, VECTOR& ipiv) {
    return LAPACKE_dgetrf(matrix_layout, m, n, storage(a), lda(a), storage(ipiv));
  }
};

template<>
struct getrf_dispatch<std::complex<float>, lapack_int> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int getrf(int matrix_layout, lapack_int m, lapack_int n, MATRIX& a, VECTOR& ipiv) {
    return LAPACKE_cgetrf(matrix_layout, m, n, storage(a), lda(a), storage(ipiv));
  }
};

template<>
struct getrf_dispatch<std::complex<double>, lapack_int> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int getrf(int matrix_layout, lapack_int m, lapack_int n, MATRIX& a, VECTOR& ipiv) {
    return LAPACKE_zgetrf(matrix_layout, m, n, storage(a), lda(a), storage(ipiv));
  }
};

}
  
template<typename MATRIX, typename VECTOR>
lapack_int getrf(MATRIX& a, VECTOR& ipiv) {
  lapack_int m = rows(a);
  lapack_int n = cols(a);
  if (size(ipiv) != std::min(m, n))
    throw std::invalid_argument("vector ipiv size mismatch");
  return getrf_dispatch<typename matrix_traits<MATRIX>::value_type,
                        typename vector_traits<VECTOR>::value_type
                        >::getrf((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR),
                                 m, n, a, ipiv);
}

} // end namespace lapackx
} // end namespace rokko

#endif // ROKKO_LAPACKX_GETRF_HPP
