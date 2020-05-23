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

#ifndef ROKKO_SCALAPACK_PSYEVR_HPP
#define ROKKO_SCALAPACK_PSYEVR_HPP

#include <rokko/cscalapack.h>
#include <rokko/lapack/storage.hpp>
#include <rokko/lapack/complex_cast.hpp>
#include <rokko/traits/real_t.hpp>
#include <rokko/traits/value_t.hpp>

namespace rokko {
namespace scalapack {

using rokko::lapack::storage;
using rokko::lapack::complex_cast;

namespace {

inline int psyevr_dispatch(char jobz, char range, char uplo, int n, float* A, int ia, int ja, const int* descA,
                           float vl, float vu, int il, int iu,
                           int& m, int& nz,
                           float* w, float* Z, int iz, int jz, const int* descZ,
                           float* work, int lwork, int* iwork, int liwork) {
  return cscalapack_pssyevr_work(jobz, range, uplo, n, A, ia, ja, descA, vl, vu, il, iu, &m, &nz, w, Z, iz, jz, descZ,
                                 work, lwork, iwork, liwork);
}

inline int psyevr_dispatch(char jobz, char range, char uplo, int n, float* A, int ia, int ja, const int* descA,
                           float vl, float vu, int il, int iu,
                           int& m, int& nz,
                           float* w, float* Z, int iz, int jz, const int* descZ) {
  return cscalapack_pssyevr(jobz, range, uplo, n, A, ia, ja, descA, vl, vu, il, iu, &m, &nz, w, Z, iz, jz, descZ);
}

inline int psyevr_dispatch(char jobz, char range, char uplo, int n, double* A, int ia, int ja, const int* descA,
                           double vl, double vu, int il, int iu,
                           int& m, int& nz,
                           double* w, double* Z, int iz, int jz, const int* descZ,
                           double* work, int lwork, int* iwork, int liwork) {
  return cscalapack_pdsyevr_work(jobz, range, uplo, n, A, ia, ja, descA, vl, vu, il, iu, &m, &nz, w, Z, iz, jz, descZ,
                                 work, lwork, iwork, liwork);
}

inline int psyevr_dispatch(char jobz, char range, char uplo, int n, double* A, int ia, int ja, const int* descA,
                           double vl, double vu, int il, int iu,
                           int& m, int& nz,
                           double* w, double* Z, int iz, int jz, const int* descZ) {
  return cscalapack_pdsyevr(jobz, range, uplo, n, A, ia, ja, descA, vl, vu, il, iu, &m, &nz, w, Z, iz, jz, descZ);
}

inline int psyevr_dispatch(char jobz, char range, char uplo, int n, std::complex<float>* A, int ia, int ja, const int* descA,
                           float vl, float vu, int il, int iu,
                           int& m, int& nz,
                           float* w, std::complex<float>* Z, int iz, int jz, const int* descZ) {
  return cscalapack_pcheevr(jobz, range, uplo, n, complex_cast(A), ia, ja, descA, vl, vu, il, iu, &m, &nz, w, complex_cast(Z), iz, jz, descZ);
}

inline int psyevr_dispatch(char jobz, char range, char uplo, int n, std::complex<double>* A, int ia, int ja, const int* descA,
                           double vl, double vu, int il, int iu,
                           int& m, int& nz,
                           double* w, std::complex<double>* Z, int iz, int jz, const int* descZ) {
  return cscalapack_pzheevr(jobz, range, uplo, n, complex_cast(A), ia, ja, descA, vl, vu, il, iu, &m, &nz, w, complex_cast(Z), iz, jz, descZ);
}

} // end of anonymous namespace

template<typename T, typename MATRIX, typename VECTOR>
int psyevr(char jobz, char range, char uplo, MATRIX& a,
           T vl, T vu, int il, int iu,
           int& m, int& nz,
           VECTOR& w, MATRIX& z) {
  static_assert(std::is_same<real_t<MATRIX>, value_t<VECTOR>>::value, "");
  static_assert(std::is_same<value_t<VECTOR>, T>::value, "");

  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  const int* descZ = z.get_mapping().get_blacs_descriptor().data();
  return psyevr_dispatch(jobz, range, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                         vl, vu, il, iu, m, nz,
                         storage(w), z.get_array_pointer(), 0, 0, descZ);
}

template<typename T, typename MATRIX, typename VECTOR0, typename VECTOR1>
int psyevr(char jobz, char range, char uplo, MATRIX& a,
           T vl, T vu, int il, int iu,
           int& m, int& nz,
           VECTOR0& w, MATRIX& z, VECTOR0& work, VECTOR1& iwork) {
  static_assert(std::is_same<real_t<MATRIX>, value_t<VECTOR0>>::value, "");
  static_assert(std::is_same<value_t<VECTOR0>, T>::value, "");

  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  const int* descZ = z.get_mapping().get_blacs_descriptor().data();
  return psyevr_dispatch(jobz, range, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                         vl, vu, il, iu, m, nz,
                         storage(w), z.get_array_pointer(), 0, 0, descZ,
                         storage(work), work.size(), storage(iwork), iwork.size());
}

// eigenvalues & eigenvectors (without jobz)
template<typename T, typename MATRIX, typename VECTOR>
int psyevr(char range, char uplo, MATRIX& a,
           T vl, T vu, int il, int iu,
           int& m, int& nz,
           VECTOR& w, MATRIX& z) {
  return psyevr('V', range, uplo, a,
                vl, vu, il, iu, m, nz,
                w, z);
}

// eigenvalues & eigenvectors, use all (without jobz, range, vl, vu, il, iu)
template<typename MATRIX, typename VECTOR>
int psyevr(char uplo, MATRIX& a,
           int& m, int& nz,
           VECTOR& w, MATRIX& z) {
  constexpr real_t<MATRIX> real_zero = 0;
  return psyevr('A', uplo, a,
                real_zero, real_zero, 0, 0, m, nz,
                w, z);
}

// eigenvalues & eigenvectors, use vl, vu (without jobz, range, il, iu)
template<typename T, typename MATRIX, typename VECTOR>
int psyevr(char uplo, MATRIX& a,
           T vl, T vu,
           int& m, int& nz,
           VECTOR& w, MATRIX& z) {
  return psyevr('V', uplo, a,
                vl, vu, 0, 0, m, nz,
                w, z);
}

// eigenvalues & eigenvectors, use il, iu (without jobz, range, vl, vu)
template<typename T, typename MATRIX, typename VECTOR>
int psyevr(char uplo, MATRIX& a,
           int il, int iu,
           int& m, int& nz,
           VECTOR& w, MATRIX& z) {
  constexpr real_t<MATRIX> real_zero = 0;
  return psyevr('I', uplo, a,
                real_zero, real_zero, il, iu, m, nz,
                w, z);
}

// only eigenvalues (without jobz)
template<typename T, typename MATRIX, typename VECTOR>
int psyevr(char range, char uplo, MATRIX& a,
           T vl, T vu, int il, int iu,
           int& m, int& nz,
           VECTOR& w) {
  static_assert(std::is_same<real_t<MATRIX>, value_t<VECTOR>>::value, "");
  static_assert(std::is_same<value_t<VECTOR>, T>::value, "");

  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  return psyevr_dispatch('N', range, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                         vl, vu, il, iu, m, nz,
                         storage(w), NULL, 0, 0, NULL);
}

// only eigenvalues, use all (without jobz, range, vl, vu, il, iu)
template<typename MATRIX, typename VECTOR>
int psyevr(char uplo, MATRIX& a,
           int& m, int& nz,
           VECTOR& w) {
  constexpr real_t<MATRIX> real_zero = 0;
  return psyevr('A', uplo, a,
                real_zero, real_zero, 0, 0, m, nz,
                w);
}

// only eigenvalues, use vl, vu (without jobz, range, il, iu)
template<typename T, typename MATRIX, typename VECTOR>
int psyevr(char uplo, MATRIX& a,
           T vl, T vu,
           int& m, int& nz,
           VECTOR& w) {
  return psyevr('V', uplo, a,
                vl, vu, 0, 0, m, nz,
                w);
}

// only eigenvalues, use il, iu (without jobz, range, vl, vu)
template<typename T, typename MATRIX, typename VECTOR>
int psyevr(char uplo, MATRIX& a,
           int il, int iu,
           int& m, int& nz,
           VECTOR& w) {
  constexpr real_t<MATRIX> real_zero = 0;
  return psyevr('I', uplo, a,
                real_zero, real_zero, il, iu, m, nz,
                w);
}

} // end namespace scalapack
} // end namespace rokko

#endif // ROKKO_SCALAPACK_PSYEVR_HPP
