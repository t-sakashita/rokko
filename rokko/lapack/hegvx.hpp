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

#ifndef ROKKO_LAPACK_HEGVX_HPP
#define ROKKO_LAPACK_HEGVX_HPP

#include <complex>
#include <lapacke.h>
#undef I
#include <rokko/traits/real_t.hpp>
#include <rokko/traits/value_t.hpp>
#include "complex_cast.hpp"
#include <rokko/lapack/storage.hpp>
#include <rokko/alias_template_function.hpp>

namespace rokko {
namespace lapack {

namespace {

template<typename T> struct hegvx_dispatch;

template<>
struct hegvx_dispatch<float> {
  template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
  static lapack_int hegvx(int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,
                          MATRIX& a, MATRIX& b,
                          float vl, float vu, lapack_int il, lapack_int iu, float abstol,
                          lapack_int& m, VECTOR& w, VECTOR_INT& ifail) {
  return LAPACKE_ssygvx(matrix_layout, itype, jobz, uplo, n, storage(a), ld(a), storage(b), ld(b), vl, vu, il, iu, abstol, &m, storage(w), storage(ifail));
  }
};

template<>
struct hegvx_dispatch<double> {
  template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
  static lapack_int hegvx(int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,
                          MATRIX& a, MATRIX& b,
                          double vl, double vu, lapack_int il, lapack_int iu, double abstol,
                          lapack_int& m, VECTOR& w, VECTOR_INT& ifail) {
    return LAPACKE_dsygvx(matrix_layout, itype, jobz, uplo, n, storage(a), ld(a), storage(b), ld(b), vl, vu, il, iu, abstol, &m, storage(w), storage(ifail));
  }
};

template<>
struct hegvx_dispatch<std::complex<float>> {
  template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
  static lapack_int hegvx(int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,
                          MATRIX& a, MATRIX& b,
                          float vl, float vu, lapack_int il, lapack_int iu, float abstol,
                          lapack_int& m, VECTOR& w, VECTOR_INT& ifail) {
    return LAPACKE_chegvx(matrix_layout, itype, jobz, uplo, n, storage(a), ld(a), storage(b), ld(b), vl, vu, il, iu, abstol, &m, storage(w), storage(ifail));
  }
};

template<>
struct hegvx_dispatch<std::complex<double>> {
  template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
  static lapack_int hegvx(int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,
                          MATRIX& a, MATRIX& b,
                          double vl, double vu, lapack_int il, lapack_int iu, double abstol,
                          lapack_int& m, VECTOR& w, VECTOR_INT& ifail) {
    return LAPACKE_zhegvx(matrix_layout, itype, jobz, uplo, n, storage(a), ld(a), storage(b), ld(b), vl, vu, il, iu, abstol, &m, storage(w), storage(ifail));
  }
};

} // end of anonymous namespace

template<typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int hegvx(lapack_int itype, char jobz, char uplo,
                 MATRIX& a, MATRIX& b,
                 double vl, double vu, lapack_int il, lapack_int iu, double abstol,
                 lapack_int& m, VECTOR& w, VECTOR_INT& ifail) {
  static_assert(std::is_same<real_t<MATRIX>, value_t<VECTOR>>::value, "");
  lapack_int n = rows(a);
  if (rows(a) != cols(a))
    throw std::invalid_argument("matrix A size mismatch");
  //if (size(w) != n)
  //  throw std::invalid_argument("vector w size mismatch");
  return hegvx_dispatch<value_t<MATRIX>>::
    hegvx((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), itype, jobz, uplo, n, a, b, w, vl, vu, il, iu, abstol, m, w, ifail);
 }

ALIAS_TEMPLATE_FUNCTION(sygvx, hegvx);

} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_HEGVX_HPP
