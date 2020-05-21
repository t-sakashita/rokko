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

#ifndef ROKKO_SCALAPACK_PSYEV_HPP
#define ROKKO_SCALAPACK_PSYEV_HPP

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

inline int psyev_dispatch(char jobz, char uplo, int n, double* A, int ia, int ja, const int* descA,
                          double* w, double* Z, int iz, int jz, const int* descZ,
                          double* work, int lwork) {
  return cscalapack_pdsyev_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ, work, lwork);
}

inline int psyev_dispatch(char jobz, char uplo, int n, float* A, int ia, int ja, const int* descA,
                          float* w, float* Z, int iz, int jz, const int* descZ) {
  return cscalapack_pssyev(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ);
}

inline int psyev_dispatch(char jobz, char uplo, int n, double* A, int ia, int ja, const int* descA,
                          double* w, double* Z, int iz, int jz, const int* descZ) {
  return cscalapack_pdsyev(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ);
}

inline int psyev_dispatch(char jobz, char uplo, int n, std::complex<float>* A, int ia, int ja, const int* descA,
                          float* w, std::complex<float>* Z, int iz, int jz, const int* descZ) {
  return cscalapack_pcheev(jobz, uplo, n, complex_cast(A), ia, ja, descA, w, complex_cast(Z), iz, jz, descZ);
}

inline int psyev_dispatch(char jobz, char uplo, int n, std::complex<double>* A, int ia, int ja, const int* descA,
                          double* w, std::complex<double>* Z, int iz, int jz, const int* descZ) {
  return cscalapack_pzheev(jobz, uplo, n, complex_cast(A), ia, ja, descA, w, complex_cast(Z), iz, jz, descZ);
}

} // end of anonymous namespace

template<typename MATRIX, typename VECTOR>
int psyev(char jobz, char uplo, MATRIX& a, VECTOR& w, MATRIX& z) {
  static_assert(std::is_same<real_t<MATRIX>, value_t<VECTOR>>::value, "");

  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  const int* descZ = z.get_mapping().get_blacs_descriptor().data();
  return psyev_dispatch(jobz, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                        storage(w), z.get_array_pointer(), 0, 0, descZ);
}

template<typename MATRIX, typename VECTOR>
int psyev(char jobz, char uplo, MATRIX& a, VECTOR& w, MATRIX& z, VECTOR& work) {
  static_assert(std::is_same<real_t<MATRIX>, value_t<VECTOR>>::value, "");

  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  const int* descZ = z.get_mapping().get_blacs_descriptor().data();
  return psyev_dispatch(jobz, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                        storage(w), z.get_array_pointer(), 0, 0, descZ, storage(work), work.size());
}

// eigenvalues & eigenvectors (without jobz)
template<typename MATRIX, typename VECTOR>
int psyev(char uplo, MATRIX& a, VECTOR& w, MATRIX& z) {
  return psyev('V', uplo, a, w, z);
}

// only eigenvalues (without jobz)
template<typename MATRIX, typename VECTOR>
int psyev(char uplo, MATRIX& a, VECTOR& w) {
  static_assert(std::is_same<real_t<MATRIX>, value_t<VECTOR>>::value, "");

  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  return psyev_dispatch('N', uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                        storage(w), NULL, 0, 0, NULL);
}

} // end namespace scalapack
} // end namespace rokko

#endif // ROKKO_SCALAPACK_PSYEV_HPP
