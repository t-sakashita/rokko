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

#ifndef ROKKO_SCALAPACK_PSYEVX_HPP
#define ROKKO_SCALAPACK_PSYEVX_HPP

#include <rokko/cscalapack.h>
#include <rokko/eigen3.hpp>
#include <rokko/lapack/storage.hpp>
#include <rokko/lapack/complex_cast.hpp>
#include <rokko/traits/real_t.hpp>
#include <rokko/traits/value_t.hpp>

namespace rokko {
namespace scalapack {

using rokko::lapack::storage;
using rokko::lapack::complex_cast;

namespace {

inline int psyevx_dispatch(char jobz, char range, char uplo, int n, float* A, int ia, int ja, const int* descA,
                           float vl, float vu, int il, int iu,
                           float abstol, int& m, int& nz,
                           float* w, float orfac, float* Z, int iz, int jz, const int* descZ,
                           float* work, int lwork, int* iwork, int liwork,
                           int* ifail, int* iclustr, float* gap) {
  return cscalapack_pssyevx_work(jobz, range, uplo, n, A, ia, ja, descA, vl, vu, il, iu, abstol, &m, &nz, w, orfac, Z, iz, jz, descZ,
                                 work, lwork, iwork, liwork,
                                 ifail, iclustr, gap);
}

inline int psyevx_dispatch(char jobz, char range, char uplo, int n, float* A, int ia, int ja, const int* descA,
                           float vl, float vu, int il, int iu,
                           float abstol, int& m, int& nz,
                           float* w, float orfac, float* Z, int iz, int jz, const int* descZ,
                           int* ifail, int* iclustr, float* gap) {
  return cscalapack_pssyevx(jobz, range, uplo, n, A, ia, ja, descA, vl, vu, il, iu, abstol, &m, &nz, w, orfac, Z, iz, jz, descZ,
                            ifail, iclustr, gap);
}

inline int psyevx_dispatch(char jobz, char range, char uplo, int n, double* A, int ia, int ja, const int* descA,
                           double vl, double vu, int il, int iu,
                           double abstol, int& m, int& nz,
                           double* w, double orfac, double* Z, int iz, int jz, const int* descZ,
                           double* work, int lwork, int* iwork, int liwork,
                           int* ifail, int* iclustr, double* gap) {
  return cscalapack_pdsyevx_work(jobz, range, uplo, n, A, ia, ja, descA, vl, vu, il, iu, abstol, &m, &nz, w, orfac, Z, iz, jz, descZ,
                                 work, lwork, iwork, liwork,
                                 ifail, iclustr, gap);
}

inline int psyevx_dispatch(char jobz, char range, char uplo, int n, double* A, int ia, int ja, const int* descA,
                           double vl, double vu, int il, int iu,
                           double abstol, int& m, int& nz,
                           double* w, double orfac, double* Z, int iz, int jz, const int* descZ,
                           int* ifail, int* iclustr, double* gap) {
  return cscalapack_pdsyevx(jobz, range, uplo, n, A, ia, ja, descA, vl, vu, il, iu, abstol, &m, &nz, w, orfac, Z, iz, jz, descZ,
                            ifail, iclustr, gap);
}

inline int psyevx_dispatch(char jobz, char range, char uplo, int n, std::complex<float>* A, int ia, int ja, const int* descA,
                           float vl, float vu, int il, int iu,
                           float abstol, int& m, int& nz,
                           float* w, float orfac, std::complex<float>* Z, int iz, int jz, const int* descZ,
                           int* ifail, int* iclustr, float* gap) {
  return cscalapack_pcheevx(jobz, range, uplo, n, complex_cast(A), ia, ja, descA, vl, vu, il, iu, abstol, &m, &nz, w, orfac, complex_cast(Z), iz, jz, descZ,
                            ifail, iclustr, gap);
}

inline int psyevx_dispatch(char jobz, char range, char uplo, int n, std::complex<double>* A, int ia, int ja, const int* descA,
                           double vl, double vu, int il, int iu,
                           double abstol, int& m, int& nz,
                           double* w, double orfac, std::complex<double>* Z, int iz, int jz, const int* descZ,
                           int* ifail, int* iclustr, double* gap) {
  return cscalapack_pzheevx(jobz, range, uplo, n, complex_cast(A), ia, ja, descA, vl, vu, il, iu, abstol, &m, &nz, w, orfac, complex_cast(Z), iz, jz, descZ,
                            ifail, iclustr, gap);
}

} // end of anonymous namespace

template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT, typename VECTOR2>
int psyevx(char jobz, char range, char uplo, MATRIX& a,
           T vl, T vu, int il, int iu,
           T abstol, int& m, int& nz,
           VECTOR& w, double orfac, MATRIX& z,
           VECTOR_INT& ifail, VECTOR_INT& iclustr, VECTOR2& gap) {
  static_assert(std::is_same_v<real_t<MATRIX>, value_t<VECTOR>>);
  static_assert(std::is_same_v<value_t<VECTOR>, value_t<VECTOR2>>);
  static_assert(std::is_same_v<value_t<VECTOR>, T>);

  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  const int* descZ = z.get_mapping().get_blacs_descriptor().data();
  return psyevx_dispatch(jobz, range, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                         vl, vu, il, iu, abstol, m, nz,
                         storage(w), orfac, z.get_array_pointer(), 0, 0, descZ,
                         storage(ifail), storage(iclustr), storage(gap));
}

template<typename T, typename MATRIX, typename VECTOR0, typename VECTOR1, typename VECTOR_INT, typename VECTOR2>
int psyevx(char jobz, char range, char uplo, MATRIX& a,
           T vl, T vu, int il, int iu,
           T abstol, int& m, int& nz,
           VECTOR0& w, double orfac, MATRIX& z, VECTOR0& work, VECTOR1& iwork,
           VECTOR_INT& ifail, VECTOR_INT& iclustr, VECTOR2& gap) {
  static_assert(std::is_same_v<real_t<MATRIX>, value_t<VECTOR0>>);
  static_assert(std::is_same_v<value_t<VECTOR0>, value_t<VECTOR2>>);
  static_assert(std::is_same_v<value_t<VECTOR0>, T>);

  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  const int* descZ = z.get_mapping().get_blacs_descriptor().data();
  return psyevx_dispatch(jobz, range, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                         vl, vu, il, iu, abstol, m, nz,
                         storage(w), orfac, z.get_array_pointer(), 0, 0, descZ,
                         storage(work), work.size(), storage(iwork), iwork.size(),
                         storage(ifail), storage(iclustr), storage(gap));
}

// eigenvalues & eigenvectors (without jobz)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT, typename VECTOR2>
int psyevx(char range, char uplo, MATRIX& a,
           T vl, T vu, int il, int iu,
           T abstol, int& m, int& nz,
           VECTOR& w, T orfac, MATRIX& z,
           VECTOR_INT& ifail, VECTOR_INT& iclustr, VECTOR2& gap) {
  return psyevx('V', range, uplo, a,
                vl, vu, il, iu, abstol, m, nz,
                w, orfac, z,
                ifail, iclustr, gap);
}

// eigenvalues & eigenvectors, use all (without jobz, range, vl, vu, il, iu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT, typename VECTOR2>
int psyevx(char uplo, MATRIX& a,
           T abstol, int& m, int& nz,
           VECTOR& w, T orfac, MATRIX& z,
           VECTOR_INT& ifail, VECTOR_INT& iclustr, VECTOR2& gap) {
  constexpr real_t<MATRIX> real_zero = 0;
  return psyevx('A', uplo, a,
                real_zero, real_zero, 0, 0, abstol, m, nz,
                w, orfac, z,
                ifail, iclustr, gap);
}

// eigenvalues & eigenvectors, use vl, vu (without jobz, range, il, iu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT, typename VECTOR2>
int psyevx(char uplo, MATRIX& a,
           T vl, T vu,
           T abstol, int& m, int& nz,
           VECTOR& w, T orfac, MATRIX& z,
           VECTOR_INT& ifail, VECTOR_INT& iclustr, VECTOR2& gap) {
  return psyevx('V', uplo, a,
                vl, vu, 0, 0, abstol, m, nz,
                w, orfac, z,
                ifail, iclustr, gap);
}

// eigenvalues & eigenvectors, use il, iu (without jobz, range, vl, vu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT, typename VECTOR2>
int psyevx(char uplo, MATRIX& a,
           int il, int iu,
           T abstol, int& m, int& nz,
           VECTOR& w, T orfac, MATRIX& z,
           VECTOR_INT& ifail, VECTOR_INT& iclustr, VECTOR2& gap) {
  constexpr real_t<MATRIX> real_zero = 0;
  return psyevx('I', uplo, a,
                real_zero, real_zero, il, iu, abstol, m, nz,
                w, orfac, z,
                ifail, iclustr, gap);
}

// only eigenvalues (without jobz)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT, typename VECTOR2>
int psyevx(char range, char uplo, MATRIX& a,
           T vl, T vu, int il, int iu,
           T abstol, int& m, int& nz,
           VECTOR& w, T orfac,
           VECTOR_INT& ifail, VECTOR_INT& iclustr, VECTOR2& gap) {
  static_assert(std::is_same_v<real_t<MATRIX>, value_t<VECTOR>>);
  static_assert(std::is_same_v<value_t<VECTOR>, value_t<VECTOR2>>);
  static_assert(std::is_same_v<value_t<VECTOR>, T>);

  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  return psyevx_dispatch('N', range, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                         vl, vu, il, iu, abstol, m, nz,
                         storage(w), orfac, NULL, 0, 0, NULL,
                         storage(ifail), storage(iclustr), storage(gap));
}

// only eigenvalues, use all (without jobz, range, vl, vu, il, iu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT, typename VECTOR2>
int psyevx(char uplo, MATRIX& a,
           T abstol, int& m, int& nz,
           VECTOR& w, T orfac,
           VECTOR_INT& ifail, VECTOR_INT& iclustr, VECTOR2& gap) {
  constexpr real_t<MATRIX> real_zero = 0;
  return psyevx('A', uplo, a,
                real_zero, real_zero, 0, 0, abstol, m, nz,
                w, orfac,
                ifail, iclustr, gap);
}

// only eigenvalues, use vl, vu (without jobz, range, il, iu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT, typename VECTOR2>
int psyevx(char uplo, MATRIX& a,
           T vl, T vu,
           T abstol, int& m, int& nz,
           VECTOR& w, T orfac,
           VECTOR_INT& ifail, VECTOR_INT& iclustr, VECTOR2& gap) {
  return psyevx('V', uplo, a,
                vl, vu, 0, 0, abstol, m, nz,
                w, orfac,
                ifail, iclustr, gap);
}

// only eigenvalues, use il, iu (without jobz, range, vl, vu)
template<typename T, typename MATRIX, typename VECTOR, typename VECTOR_INT, typename VECTOR2>
int psyevx(char uplo, MATRIX& a,
           int il, int iu,
           T abstol, int& m, int& nz,
           VECTOR& w, T orfac,
           VECTOR_INT& ifail, VECTOR_INT& iclustr, VECTOR2& gap) {
  constexpr real_t<MATRIX> real_zero = 0;
  return psyevx('I', uplo, a,
                real_zero, real_zero, il, iu, abstol, m, nz,
                w, orfac,
                ifail, iclustr, gap);
}

} // end namespace scalapack
} // end namespace rokko

#endif // ROKKO_SCALAPACK_PSYEVX_HPP
