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

#ifndef ROKKO_SCALAPACK_PSYEVD_HPP
#define ROKKO_SCALAPACK_PSYEVD_HPP

#include <rokko/cscalapack.h>
#include <rokko/lapack/complex_cast.hpp>
#include <rokko/traits/real_t.hpp>
#include <rokko/traits/value_t.hpp>

namespace rokko {
namespace scalapack {

using rokko::lapack::storage;
using rokko::lapack::complex_cast;

namespace {

inline int psyevd_dispatch(char jobz, char uplo, int n, float* A, int ia, int ja, const int* descA,
                           float* w, float* Z, int iz, int jz, const int* descZ,
                           float* work, int lwork, int* iwork, int liwork) {
  return cscalapack_pssyevd_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ,
                                 work, lwork, iwork, liwork);
}

inline int psyevd_dispatch(char jobz, char uplo, int n, float* A, int ia, int ja, const int* descA,
                           float* w, float* Z, int iz, int jz, const int* descZ) {
  return cscalapack_pssyevd(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ);
}

inline int psyevd_dispatch(char jobz, char uplo, int n, double* A, int ia, int ja, const int* descA,
                           double* w, double* Z, int iz, int jz, const int* descZ,
                           double* work, int lwork, int* iwork, int liwork) {
  return cscalapack_pdsyevd_work(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ,
                                 work, lwork, iwork, liwork);
}

inline int psyevd_dispatch(char jobz, char uplo, int n, double* A, int ia, int ja, const int* descA,
                           double* w, double* Z, int iz, int jz, const int* descZ) {
  return cscalapack_pdsyevd(jobz, uplo, n, A, ia, ja, descA, w, Z, iz, jz, descZ);
}

inline int psyevd_dispatch(char jobz, char uplo, int n, std::complex<float>* A, int ia, int ja, const int* descA,
                           float* w, std::complex<float>* Z, int iz, int jz, const int* descZ) {
  return cscalapack_pcheevd(jobz, uplo, n, complex_cast(A), ia, ja, descA, w, complex_cast(Z), iz, jz, descZ);
}

inline int psyevd_dispatch(char jobz, char uplo, int n, std::complex<double>* A, int ia, int ja, const int* descA,
                           double* w, std::complex<double>* Z, int iz, int jz, const int* descZ) {
  return cscalapack_pzheevd(jobz, uplo, n, complex_cast(A), ia, ja, descA, w, complex_cast(Z), iz, jz, descZ);
}

} // end of anonymous namespace

template<typename MATRIX, typename VECTOR>
int psyevd(char jobz, char uplo, MATRIX& a, VECTOR& w, MATRIX& z) {
  static_assert(std::is_same<real_t<MATRIX>, value_t<VECTOR>>::value);

  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  const int* descZ = z.get_mapping().get_blacs_descriptor().data();
  return psyevd_dispatch(jobz, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                        storage(w), z.get_array_pointer(), 0, 0, descZ);
}

template<typename MATRIX, typename VECTOR0, typename VECTOR1>
int psyevd(char jobz, char uplo, MATRIX& a, VECTOR0& w, MATRIX& z, VECTOR0& work, VECTOR1& iwork) {
  static_assert(std::is_same<real_t<MATRIX>, value_t<VECTOR0>>::value);

  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  const int* descZ = z.get_mapping().get_blacs_descriptor().data();
  return psyevd_dispatch(jobz, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                         storage(w), z.get_array_pointer(), 0, 0, descZ,
                         storage(work), work.size(), storage(iwork), iwork.size());
}

// eigenvalues & eigenvectors (without jobz)
template<typename MATRIX, typename VECTOR>
int psyevd(char uplo, MATRIX& a, VECTOR& w, MATRIX& z) {
  return psyevd('V', uplo, a, w, z);
}

} // end namespace scalapack
} // end namespace rokko

#endif // ROKKO_SCALAPACK_PSYEVD_HPP
