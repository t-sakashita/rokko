/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2017-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_LAPACK_HTEIN_HPP
#define ROKKO_LAPACK_HTEIN_HPP

#include <complex>
#include <lapacke.h>
#undef I
#include <rokko/traits/real_t.hpp>
#include <rokko/traits/value_t.hpp>
#include "complex_cast.hpp"
#include <rokko/alias_template_function.hpp>

namespace rokko {
namespace lapack {

namespace {

template<typename T> struct htein_dispatch;

template<>
struct htein_dispatch<float> {
  template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
  static lapack_int htein(int matrix_layout, lapack_int n,
                          VECTOR& d, VECTOR& e,
                          lapack_int m, VECTOR& w,
                          MATRIX& z,
                          const VECTOR_INT& iblock, const VECTOR_INT& isplit,
                          VECTOR_INT& ifailv) {
    return LAPACKE_sstein(matrix_layout, n, storage(d), storage(e),
                          m, storage(w), storage(iblock), storage(isplit), storage(z), ld(z), storage(ifailv));
  }
};

template<>
struct htein_dispatch<double> {
  template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
  static lapack_int htein(int matrix_layout, lapack_int n,
                          VECTOR& d, VECTOR& e,
                          lapack_int m, VECTOR& w,
                          MATRIX& z,
                          const VECTOR_INT& iblock, const VECTOR_INT& isplit,
                          VECTOR_INT& ifailv) {
    return LAPACKE_dstein(matrix_layout, n, storage(d), storage(e),
                          m, storage(w), storage(iblock), storage(isplit), storage(z), ld(z), storage(ifailv));
  }
};

template<>
struct htein_dispatch<std::complex<float>> {
  template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
  static lapack_int htein(int matrix_layout, lapack_int n,
                          VECTOR& d, VECTOR& e,
                          lapack_int m, VECTOR& w,
                          MATRIX& z,
                          const VECTOR_INT& iblock, const VECTOR_INT& isplit,
                          VECTOR_INT& ifailv) {
    return LAPACKE_cstein(matrix_layout, n, storage(d), storage(e),
                          m, storage(w), srorage(iblock), storage(isplit), complex_cast(storage(z)), ld(z), storage(ifailv));
  }
};

template<>
struct htein_dispatch<std::complex<double>> {
  template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
  static lapack_int htein(int matrix_layout, lapack_int n,
                          VECTOR& d, VECTOR& e,
                          lapack_int m, VECTOR& w,
                          MATRIX& z,
                          const VECTOR_INT& iblock, const VECTOR_INT& isplit,
                          VECTOR_INT& ifailv) {
    return LAPACKE_zstein(matrix_layout, n, storage(d), storage(e),
                          m, storage(w), storage(iblock), storage(isplit), complex_cast(storage(z)), ld(z), storage(ifailv));
  }
};

} // end of anonymous namespace

template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int htein(VECTOR& d, VECTOR& e, lapack_int m, VECTOR& w, MATRIX& z,
                 const VECTOR_INT& iblock, const VECTOR_INT& isplit, VECTOR_INT& ifailv) {
  static_assert(std::is_same<real_t<MATRIX>, value_t<VECTOR>>::value, "");
  lapack_int n = size(d);
  if (rows(z) != n)
    throw std::invalid_argument("matrix Z row size mismatch");
  if (cols(z) < m)
    throw std::invalid_argument("matrix Z col size mismatch");
  if (size(d) != n)
    throw std::invalid_argument("vector d size mismatch");
  if (size(e) != (n-1))
    throw std::invalid_argument("vector e size mismatch");
  if (size(w) != n)
    throw std::invalid_argument("vector w size mismatch");
  return htein_dispatch<value_t<MATRIX>>::
    htein((is_col_major(z) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), n, d, e, m, w, z, iblock, isplit, ifailv);
}

ALIAS_TEMPLATE_FUNCTION(stein, htein);

} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_HTEIN_HPP
