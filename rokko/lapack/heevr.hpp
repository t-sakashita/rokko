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

#pragma once

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

template<typename T> struct heevr_dispatch;

template<>
struct heevr_dispatch<float> {
  template<typename MATRIX, typename MATRIX2, typename VECTOR, typename VECTOR_INT>
  static lapack_int heevr(int matrix_layout, char jobz, char range, char uplo, lapack_int n, MATRIX& a,
                          float vl, float vu, lapack_int il, lapack_int iu, float abstol,
                          lapack_int& m, VECTOR& w, MATRIX2& z, VECTOR_INT& isuppz) {
    return LAPACKE_ssyevr(matrix_layout, jobz, range, uplo, n, storage(a), ld(a), vl, vu, il, iu, abstol, &m, storage(w), storage(z), ld(z), storage(isuppz));
  }
};

template<>
struct heevr_dispatch<double> {
  template<typename MATRIX, typename MATRIX2, typename VECTOR, typename VECTOR_INT>
  static lapack_int heevr(int matrix_layout, char jobz, char range, char uplo, lapack_int n, MATRIX& a,
                          double vl, double vu, lapack_int il, lapack_int iu, double abstol,
                          lapack_int& m, VECTOR& w, MATRIX2& z, VECTOR_INT& isuppz) {
    return LAPACKE_dsyevr(matrix_layout, jobz, range, uplo, n, storage(a), ld(a), vl, vu, il, iu, abstol, &m, storage(w), storage(z), ld(z), storage(isuppz));
  }
};

template<>
struct heevr_dispatch<std::complex<float>> {
  template<typename MATRIX, typename MATRIX2, typename VECTOR, typename VECTOR_INT>
  static lapack_int heevr(int matrix_layout, char jobz, char range, char uplo, lapack_int n, MATRIX& a,
                          float vl, float vu, lapack_int il, lapack_int iu, float abstol,
                          lapack_int& m, VECTOR& w, MATRIX2& z, VECTOR_INT& isuppz) {
    return LAPACKE_cheevr(matrix_layout, jobz, range, uplo, n, complex_cast(storage(a)), ld(a), vl, vu, il, iu, abstol, &m, storage(w), complex_cast(storage(z)), ld(z), storage(isuppz));
  }
};

template<>
struct heevr_dispatch<std::complex<double>> {
  template<typename MATRIX, typename MATRIX2, typename VECTOR, typename VECTOR_INT>
  static lapack_int heevr(int matrix_layout, char jobz, char range, char uplo, lapack_int n, MATRIX& a,
                          double vl, double vu, lapack_int il, lapack_int iu, double abstol,
                          lapack_int& m, VECTOR& w, MATRIX2& z, VECTOR_INT& isuppz) {
    return LAPACKE_zheevr(matrix_layout, jobz, range, uplo, n, complex_cast(storage(a)), ld(a), vl, vu, il, iu, abstol, &m, storage(w), complex_cast(storage(z)), ld(z), storage(isuppz));
  }
};

} // end of anonymous namespace

template<typename T, typename MATRIX, typename MATRIX2, typename VECTOR, typename VECTOR_INT>
lapack_int heevr(char jobz, char range, char uplo, MATRIX& a,
                 T vl, T vu, lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, MATRIX2& z, VECTOR_INT& isuppz) {
  static_assert(std::is_same_v<real_t<MATRIX>, value_t<VECTOR>>);
  lapack_int n = rows(a);
  if (rows(a) != cols(a))
    throw std::invalid_argument("matrix A size mismatch");
  //if (size(w) != n)
  //  throw std::invalid_argument("vector w size mismatch");
  return heevr_dispatch<value_t<MATRIX>>::
    heevr((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), jobz, range, uplo, n, a, vl, vu, il, iu, abstol, m, w, z, isuppz);
}

// eigenvalues & eigenvectors (without jobz)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevr(char range, char uplo, MATRIX& a,
                 T vl, T vu, lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& isuppz) {
  return heevr('V', range, uplo, a, vl, vu, il, iu, abstol, m, w, z, isuppz);
}

// eigenvalues & eigenvectors, use all (without jobz, range, vl, vu, il, iu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevr(char uplo, MATRIX& a,
                 T abstol,
                 lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& isuppz) {
  constexpr real_t<MATRIX> real_zero = 0;
  return heevr('A', uplo, a, real_zero, real_zero, 0, 0, abstol, m, w, z, isuppz);
}

// eigenvalues & eigenvectors, use vl, vu (without jobz, range, il, iu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevr(char uplo, MATRIX& a,
                 T vl, T vu, T abstol,
                 lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& isuppz) {
  return heevr('V', uplo, a, vl, vu, 0, 0, abstol, m, w, z, isuppz);
}

// eigenvalues & eigenvectors, use il, iu (without jobz, range, vl, vu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevr(char uplo, MATRIX& a,
                 lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, MATRIX& z, VECTOR_INT& isuppz) {
  constexpr real_t<MATRIX> real_zero = 0;
  return heevr('I', uplo, a, real_zero, real_zero, il, iu, abstol, m, w, z, isuppz);
}

// only eigenvalues (without jobz)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevr(char range, char uplo, MATRIX& a,
                 T vl, T vu, lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, VECTOR_INT& isuppz) {
  constexpr null_matrix<value_t<MATRIX>> z_null;
  return heevr('N', range, uplo, a, vl, vu, il, iu, abstol, m, w, z_null, isuppz);
}

// only eigenvalues, use all (without jobz, range, vl, vu, il, iu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevr(char uplo, MATRIX& a,
                 T abstol,
                 lapack_int& m, VECTOR& w, VECTOR_INT& isuppz) {
  constexpr real_t<MATRIX> real_zero = 0;
  return heevr('A', uplo, a, real_zero, real_zero, 0, 0, abstol, m, w, isuppz);
}

// only eigenvalues, use vl, vu (without jobz, range, il, iu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevr(char uplo, MATRIX& a,
                 T vl, T vu, T abstol,
                 lapack_int& m, VECTOR& w, VECTOR_INT& isuppz) {
  return heevr('V', uplo, a, vl, vu, 0, 0, abstol, m, w, isuppz);
}

// only eigenvalues, use il, iu (without jobz, range, vl, vu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT>
lapack_int heevr(char uplo, MATRIX& a,
                 lapack_int il, lapack_int iu, T abstol,
                 lapack_int& m, VECTOR& w, VECTOR_INT& isuppz) {
  constexpr real_t<MATRIX> real_zero = 0;
  return heevr('I', uplo, a, real_zero, real_zero, il, iu, abstol, m, w, isuppz);
}

ALIAS_TEMPLATE_FUNCTION(syevr, heevr);

} // end namespace lapack
} // end namespace rokko
