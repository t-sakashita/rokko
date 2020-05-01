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

#ifndef ROKKO_LAPACK_HEEV_HPP
#define ROKKO_LAPACK_HEEV_HPP

#include <complex>
#include <lapacke.h>
#undef I
#include <rokko/traits/norm_t.hpp>
#include <rokko/traits/value_t.hpp>
#include "complex_cast.hpp"

namespace rokko {
namespace lapack {

namespace {

template<typename T> struct heev_dispatch;
  
template<>
struct heev_dispatch<float> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int heev(int matrix_layout, char jobz, char uplo, lapack_int n,
                         MATRIX& a, VECTOR& w) {
    return LAPACKE_ssyev(matrix_layout, jobz, uplo, n, storage(a), ld(a), storage(w));
  }
};
  
template<>
struct heev_dispatch<double> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int heev(int matrix_layout, char jobz, char uplo, lapack_int n,
                         MATRIX& a, VECTOR& w) {
    return LAPACKE_dsyev(matrix_layout, jobz, uplo, n, storage(a), ld(a), storage(w));
  }
};
  
template<>
struct heev_dispatch<std::complex<float>> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int heev(int matrix_layout, char jobz, char uplo, lapack_int n,
                         MATRIX& a, VECTOR& w) {
    return LAPACKE_cheev(matrix_layout, jobz, uplo, n, complex_cast(storage(a)), ld(a),
                         storage(w));
  }
};
  
template<>
struct heev_dispatch<std::complex<double>> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int heev(int matrix_layout, char jobz, char uplo, lapack_int n,
                         MATRIX& a, VECTOR& w) {
    return LAPACKE_zheev(matrix_layout, jobz, uplo, n, complex_cast(storage(a)), ld(a),
                         storage(w));
  }
};

}
  
template<typename MATRIX, typename VECTOR>
lapack_int heev(char jobz, char uplo, MATRIX& a, VECTOR& w) {
  BOOST_STATIC_ASSERT(std::is_same<norm_t<MATRIX>, value_t<VECTOR>>::value);
  lapack_int n = rows(a);
  if (rows(a) != cols(a))
    throw std::invalid_argument("matrix A size mismatch");
  if (size(w) != n)
    throw std::invalid_argument("vector w size mismatch");
  return heev_dispatch<value_t<MATRIX>>::
    heev((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), jobz, uplo, n, a, w);
}
  
template<typename MATRIX, typename VECTOR>
lapack_int syev(char jobz, char uplo, MATRIX& a, VECTOR& w) {
  return heev(jobz, uplo, a, w);
};

template<typename MATRIX, typename VECTOR0, typename VECTOR1>
lapack_int syev(char jobz, char uplo, MATRIX& a, VECTOR0& w, VECTOR1& work,
                VECTOR1& rwork) {
  return heev(jobz, uplo, a, w, work, rwork);
};

} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_HEEV_HPP
