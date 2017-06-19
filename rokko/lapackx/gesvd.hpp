/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2017 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_LAPACKX_GESVD_HPP
#define ROKKO_LAPACKX_GESVD_HPP

#include <complex>
#include <stdexcept>
#include <lapacke.h>

namespace rokko {
namespace lapackx {

namespace {

template<typename T1, typename T2> struct gesvd_dispatch;
  
template<>
struct gesvd_dispatch<float, float> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int gesvd(int matrix_layout, char jobu, char jobvt,
                          lapack_int m, lapack_int n, MATRIX& a, VECTOR& s,
                          MATRIX& u, MATRIX& vt, VECTOR& work) {
    return LAPACKE_sgesvd_work(matrix_layout, jobu, jobvt, m, n, 
                               storage(a), lda(a), storage(s), storage(u), lda(u),
                               storage(vt), lda(vt), storage(work), size(work));
  }
};
  
template<>
struct gesvd_dispatch<double, double> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int gesvd(int matrix_layout, char jobu, char jobvt,
                          lapack_int m, lapack_int n, MATRIX& a, VECTOR& s,
                          MATRIX& u, MATRIX& vt, VECTOR& work) {
    return LAPACKE_dgesvd_work(matrix_layout, jobu, jobvt, m, n, 
                               storage(a), lda(a), storage(s), storage(u), lda(u),
                               storage(vt), lda(vt), storage(work), size(work));
  }
};
  
template<>
struct gesvd_dispatch<std::complex<float>, std::complex<float> > {
  template<typename MATRIX, typename VECTOR>
  static lapack_int gesvd(int matrix_layout, char jobu, char jobvt,
                          lapack_int m, lapack_int n, MATRIX& a, VECTOR& s,
                          MATRIX& u, MATRIX& vt, VECTOR& work) {
    return LAPACKE_cgesvd_work(matrix_layout, jobu, jobvt, m, n, 
                               storage(a), lda(a), storage(s), storage(u), lda(u),
                               storage(vt), lda(vt), storage(work), size(work));
  }
};
  
template<>
struct gesvd_dispatch<std::complex<double>, std::complex<double> > {
  template<typename MATRIX, typename VECTOR>
  static lapack_int gesvd(int matrix_layout, char jobu, char jobvt,
                          lapack_int m, lapack_int n, MATRIX& a, VECTOR& s,
                          MATRIX& u, MATRIX& vt, VECTOR& work) {
    return LAPACKE_zgesvd_work(matrix_layout, jobu, jobvt, m, n, 
                               storage(a), lda(a), storage(s), storage(u), lda(u),
                               storage(vt), lda(vt), storage(work), size(work));
  }
};

}
  
// FIXME: gesvd

template<typename MATRIX, typename VECTOR>
lapack_int gesvd_work(char jobu, char jobvt, MATRIX& a, VECTOR& s,
                      MATRIX& u, MATRIX& vt, VECTOR& work) {
  lapack_int m = rows(a);
  lapack_int n = cols(a);
  lapack_int r = std::min(m, n);
  if (size(s) != r)
    throw std::invalid_argument("vector S size mismatch");
  if (jobu == 'A' && (rows(u) != m || cols(u) != m))
    throw std::invalid_argument("matrix U size mismatch");
  if (jobu == 'S' && (rows(u) != m || cols(u) != r))
    throw std::invalid_argument("matrix U size mismatch");
  if (jobvt == 'A' && (rows(vt) != n || cols(vt) != n))
    throw std::invalid_argument("matrix Vt size mismatch");
  if (jobvt == 'S' && (rows(vt) != r || cols(vt) != n))
    throw std::invalid_argument("matrix Vt size mismatch");
  return gesvd_dispatch<typename value_t<MATRIX>::type, typename value_t<VECTOR>::type>
    ::gesvd((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR),
            jobu, jobvt, m, n, a, s, u, vt, work);
}

} // end namespace lapackx
} // end namespace rokko

#endif // ROKKO_LAPACKX_GESVD_HPP
