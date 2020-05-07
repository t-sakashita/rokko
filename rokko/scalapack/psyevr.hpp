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

namespace rokko {
namespace scalapack {

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

} // end of anonymous namespace

template<typename T, typename MATRIX, typename VECTOR>
int psyevr(char jobz, char range, char uplo, MATRIX& a,
           T vl, T vu, int il, int iu,
           int& m, int& nz,
           VECTOR& w, MATRIX& z) {
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
  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  const int* descZ = z.get_mapping().get_blacs_descriptor().data();
  return psyevr_dispatch('V', range, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                         vl, vu, il, iu, m, nz,
                         storage(w), z.get_array_pointer(), 0, 0, descZ);
}

// only eigenvalues (without jobz)
template<typename T, typename MATRIX, typename VECTOR>
int psyevr(char range, char uplo, MATRIX& a,
           T vl, T vu, int il, int iu,
           int& m, int& nz,
           VECTOR& w) {
  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  return psyevr_dispatch('N', range, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                         vl, vu, il, iu, m, nz,
                         storage(w), NULL, 0, 0, NULL);
}

} // end namespace scalapack
} // end namespace rokko

#endif // ROKKO_SCALAPACK_PSYEVR_HPP
