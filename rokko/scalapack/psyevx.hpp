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

namespace rokko {
namespace scalapack {

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

} // end of anonymous namespace

template<typename T, typename MATRIX, typename VECTOR>
int psyevx(char jobz, char range, char uplo, MATRIX& a,
           T vl, T vu, int il, int iu,
           T abstol, int& m, int& nz,
           VECTOR& w, double orfac, MATRIX& z,
           int* ifail, int* iclustr, double* gap) {
  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  const int* descZ = z.get_mapping().get_blacs_descriptor().data();
  return psyevx_dispatch(jobz, range, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                         vl, vu, il, iu, abstol, m, nz,
                         storage(w), orfac, z.get_array_pointer(), 0, 0, descZ,
                         ifail, iclustr, gap);
}

template<typename T, typename MATRIX, typename VECTOR0, typename VECTOR1>
int psyevx(char jobz, char range, char uplo, MATRIX& a,
           T vl, T vu, int il, int iu,
           T abstol, int& m, int& nz,
           VECTOR0& w, double orfac, MATRIX& z, VECTOR0& work, VECTOR1& iwork,
           int* ifail, int* iclustr, double* gap) {
  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  const int* descZ = z.get_mapping().get_blacs_descriptor().data();
  return psyevx_dispatch(jobz, range, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                         vl, vu, il, iu, abstol, m, nz,
                         storage(w), orfac, z.get_array_pointer(), 0, 0, descZ,
                         storage(work), work.size(), storage(iwork), iwork.size(),
                         ifail, iclustr, gap);
}

// eigenvalues & eigenvectors (without jobz)
template<typename T, typename MATRIX, typename VECTOR>
int psyevx(char range, char uplo, MATRIX& a,
           T vl, T vu, int il, int iu,
           T abstol, int& m, int& nz,
           VECTOR& w, T orfac, MATRIX& z,
           int* ifail, int* iclustr, T* gap) {
  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  const int* descZ = z.get_mapping().get_blacs_descriptor().data();
  return psyevx_dispatch('V', range, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                         vl, vu, il, iu, abstol, m, nz,
                         storage(w), orfac, z.get_array_pointer(), 0, 0, descZ,
                         ifail, iclustr, gap);
}

// only eigenvalues (without jobz)
template<typename T, typename MATRIX, typename VECTOR>
int psyevx(char range, char uplo, MATRIX& a,
           T vl, T vu, int il, int iu,
           T abstol, int& m, int& nz,
           VECTOR& w, T orfac,
           int* ifail, int* iclustr, T* gap) {
  const int* descA = a.get_mapping().get_blacs_descriptor().data();
  return psyevx_dispatch('N', range, uplo, a.get_m_global(), a.get_array_pointer(), 0, 0, descA,
                         vl, vu, il, iu, abstol, m, nz,
                         storage(w), orfac, NULL, 0, 0, NULL,
                         ifail, iclustr, gap);
}

} // end namespace scalapack
} // end namespace rokko

#endif // ROKKO_SCALAPACK_PSYEVX_HPP
