/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_LAPACK_HEEVX_HPP
#define ROKKO_LAPACK_HEEVX_HPP

#include <complex>
#include <lapacke.h>
#undef I
#include <rokko/traits/norm_t.hpp>
#include <rokko/traits/value_t.hpp>
#include "complex_cast.hpp"
#include <rokko/lapack/storage.hpp>
#include <rokko/alias_template_function.hpp>

namespace rokko {
namespace lapack {

namespace {

template<typename T> struct heevx_dispatch;

template<>
struct heevx_dispatch<float> {
  template<typename MATRIX, typename MATRIX2, typename VECTOR, typename VECTOR_INT>
  static lapack_int heevx(int matrix_layout, char jobz, char range, char uplo, lapack_int n, MATRIX& a,
                          float vl, float vu, lapack_int il, lapack_int iu, float abstol,
                          lapack_int& m, VECTOR& w, MATRIX2& z, VECTOR_INT& ifail) {
    return LAPACKE_ssyevx(matrix_layout, jobz, range, uplo, n, storage(a), ld(a), vl, vu, il, iu, abstol, &m, storage(w), storage(z), ld(z), storage(ifail));
  }
};

template<>
struct heevx_dispatch<double> {
  template<typename MATRIX, typename MATRIX2, typename VECTOR, typename VECTOR_INT>
  static lapack_int heevx(int matrix_layout, char jobz, char range, char uplo, lapack_int n, MATRIX& a,
                          double vl, double vu, lapack_int il, lapack_int iu, double abstol,
                          lapack_int& m, VECTOR& w, MATRIX2& z, VECTOR_INT& ifail) {
    return LAPACKE_dsyevx(matrix_layout, jobz, range, uplo, n, storage(a), ld(a), vl, vu, il, iu, abstol, &m, storage(w), storage(z), ld(z), storage(ifail));
  }
};

template<>
struct heevx_dispatch<std::complex<float>> {
  template<typename MATRIX, typename MATRIX2, typename VECTOR, typename VECTOR_INT>
  static lapack_int heevx(int matrix_layout, char jobz, char range, char uplo, lapack_int n, MATRIX& a,
                          float vl, float vu, lapack_int il, lapack_int iu, float abstol,
                          lapack_int& m, VECTOR& w, MATRIX2& z, VECTOR_INT& ifail) {
    return LAPACKE_cheevx(matrix_layout, jobz, range, uplo, n, complex_cast(storage(a)), ld(a), vl, vu, il, iu, abstol, &m, storage(w), complex_cast(storage(z)), ld(z), storage(ifail));
  }
};

template<>
struct heevx_dispatch<std::complex<double>> {
  template<typename MATRIX, typename MATRIX2, typename VECTOR, typename VECTOR_INT>
  static lapack_int heevx(int matrix_layout, char jobz, char range, char uplo, lapack_int n, MATRIX& a,
                          double vl, double vu, lapack_int il, lapack_int iu, double abstol,
                          lapack_int& m, VECTOR& w, MATRIX2& z, VECTOR_INT& ifail) {
    return LAPACKE_zheevx(matrix_layout, jobz, range, uplo, n, complex_cast(storage(a)), ld(a), vl, vu, il, iu, abstol, &m, storage(w), complex_cast(storage(z)), ld(z), storage(ifail));
  }
};

} // end of anonymous namespace

template<typename T, typename MATRIX, typename MATRIX2, typename VECTOR, typename VECTOR_INT>
lapack_int heevx(char jobz, char range, char uplo, MATRIX& a,
                 T vl, T vu, lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, MATRIX2& z, VECTOR_INT& ifail) {
  BOOST_STATIC_ASSERT(std::is_same<norm_t<MATRIX>, value_t<VECTOR>>::value);
  lapack_int n = rows(a);
  if (rows(a) != cols(a))
    throw std::invalid_argument("matrix A size mismatch");
  //if (size(w) != n)
  //  throw std::invalid_argument("vector w size mismatch");
  return heevx_dispatch<value_t<MATRIX>>::
    heevx((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), jobz, range, uplo, n, a, vl, vu, il, iu, abstol, m, w, z, ifail);
}

// eigenvalues & eigenvectors (without jobz)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevx(char range, char uplo, MATRIX& a,
                 T vl, T vu, lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& ifail) {
  return heevx('V', range, uplo, a, vl, vu, il, iu, abstol, m, w, z, ifail);
}

// eigenvalues & eigenvectors, use all (without jobz, range, vl, vu, il, iu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevx(char uplo, MATRIX& a,
                 T abstol,
                 lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& ifail) {
  return heevx('A', uplo, a, 0., 0., 0, 0, abstol, m, w, z, ifail);
}

// eigenvalues & eigenvectors, use vl, vu (without jobz, range, il, iu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevx(char uplo, MATRIX& a,
                 T vl, T vu, T abstol,
                 lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& ifail) {
  return heevx('V', uplo, a, vl, vu, 0, 0, abstol, m, w, z, ifail);
}

// eigenvalues & eigenvectors, use il, iu (without jobz, range, vl, vu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevx(char uplo, MATRIX& a,
                 lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& ifail) {
  return heevx('I', uplo, a, 0., 0., il, iu, abstol, m, w, z, ifail);
}

// only eigenvalues (without jobz)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevx(char range, char uplo, MATRIX& a,
                 T vl, T vu, lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, VECTOR_INT& ifail) {
  constexpr null_matrix<value_t<MATRIX>> z_null;
  return heevx('N', range, uplo, a, vl, vu, il, iu, abstol, m, w, z_null, ifail);
}

// only eigenvalues, use all (without jobz, range, vl, vu, il, iu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevx(char uplo, MATRIX& a,
                 T abstol,
                 lapack_int& m, VECTOR& w, VECTOR_INT& ifail) {
  return heevx('A', uplo, a, 0., 0., 0, 0, abstol, m, w, ifail);
}

// only eigenvalues, use vl, vu (without jobz, range, il, iu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevx(char uplo, MATRIX& a,
                 T vl, T vu, T abstol,
                 lapack_int& m, VECTOR& w, VECTOR_INT& ifail) {
  return heevx('V', uplo, a, vl, vu, 0, 0, abstol, m, w, ifail);
}

// only eigenvalues, use il, iu (without jobz, range, vl, vu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevx(char uplo, MATRIX& a,
                 lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, VECTOR_INT& ifail) {
  return heevx('I', uplo, a, 0., 0., il, iu, abstol, m, w, ifail);
}

ALIAS_TEMPLATE_FUNCTION(syevx, heevx);

} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_HEEVX_HPP
