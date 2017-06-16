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

#ifndef ROKKO_LAPACKX_GETRS_HPP
#define ROKKO_LAPACKX_GETRS_HPP

#include <rokko/vector_traits.hpp>
#include <rokko/matrix_traits.hpp>
#include <complex>
#include <stdexcept>
#include <lapacke.h>

namespace rokko {
namespace lapackx {

namespace {

template<typename T1, typename T2, typename T3> struct getrs_dispatch;
  
template<>
struct getrs_dispatch<float, float, lapack_int> {
  template<typename MATRIX0, typename MATRIX1, typename VECTOR>
  static lapack_int getrs(int matrix_layout, char trans, lapack_int n, lapack_int nrhs,
                          MATRIX0& a, VECTOR& ipiv, MATRIX1& b) {
    return LAPACKE_sgetrs(matrix_layout, trans, n, nrhs, storage(a), lda(a), storage(ipiv),
                          storage(b), lda(b));
  }
};

template<>
struct getrs_dispatch<double, double, lapack_int> {
  template<typename MATRIX0, typename MATRIX1, typename VECTOR>
  static lapack_int getrs(int matrix_layout, char trans, lapack_int n, lapack_int nrhs,
                          MATRIX0& a, VECTOR& ipiv, MATRIX1& b) {
    return LAPACKE_dgetrs(matrix_layout, trans, n, nrhs, storage(a), lda(a), storage(ipiv),
                          storage(b), lda(b));
  }
};

template<>
struct getrs_dispatch<std::complex<float>, std::complex<float>, lapack_int> {
  template<typename MATRIX0, typename MATRIX1, typename VECTOR>
  static lapack_int getrs(int matrix_layout, char trans, lapack_int n, lapack_int nrhs,
                          MATRIX0& a, VECTOR& ipiv, MATRIX1& b) {
    return LAPACKE_cgetrs(matrix_layout, trans, n, nrhs, storage(a), lda(a), storage(ipiv),
                          storage(b), lda(b));
  }
};

template<>
struct getrs_dispatch<std::complex<double>, std::complex<double>, lapack_int> {
  template<typename MATRIX0, typename MATRIX1, typename VECTOR>
  static lapack_int getrs(int matrix_layout, char trans, lapack_int n, lapack_int nrhs,
                          MATRIX0& a, VECTOR& ipiv, MATRIX1& b) {
    return LAPACKE_zgetrs(matrix_layout, trans, n, nrhs, storage(a), lda(a), storage(ipiv),
                          storage(b), lda(b));
  }
};

}
  
template<typename MATRIX0, typename MATRIX1, typename VECTOR>
lapack_int getrs(char trans, lapack_int nrhs, MATRIX0 const& a,
                 VECTOR const& ipiv, MATRIX1& b) {
  lapack_int n = rows(a);
  if (rows(a) != cols(a))
    throw std::invalid_argument("matrix A size mismatch");
  if (size(ipiv) != n)
    throw std::invalid_argument("vector ipiv size mismatch");
  if (rows(b) != n || cols(b) != nrhs)
    throw std::invalid_argument("matrix B size mismatch");
  return getrs_dispatch<typename matrix_traits<MATRIX0>::value_type,
                        typename matrix_traits<MATRIX1>::value_type,
                        typename vector_traits<VECTOR>::value_type
                        >::getrs((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR),
                                 trans, n, nrhs, a, ipiv, b);
}

} // end namespace lapackx
} // end namespace rokko

#endif // ROKKO_LAPACKX_GETRS_HPP
