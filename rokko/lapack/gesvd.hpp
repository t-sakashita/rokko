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

#ifndef ROKKO_LAPACK_GESVD_HPP
#define ROKKO_LAPACK_GESVD_HPP

#include <complex>
#include <stdexcept>
#include <lapacke.h>
#undef I
#include <rokko/traits/norm_t.hpp>
#include <rokko/traits/value_t.hpp>
#include "complex_cast.hpp"

namespace rokko {
namespace lapack {

namespace {

template<typename T> struct gesvd_dispatch;
  
template<>
struct gesvd_dispatch<float> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int gesvd(int matrix_layout, char jobu, char jobvt,
                          lapack_int m, lapack_int n, MATRIX& a, VECTOR& s,
                          MATRIX& u, MATRIX& vt, VECTOR& work) {
    if (size(work) == std::min(m, n) - 1)
      return LAPACKE_sgesvd(matrix_layout, jobu, jobvt, m, n, 
                            storage(a), lda(a), storage(s), storage(u), lda(u),
                            storage(vt), lda(vt), storage(work));
    else
      return LAPACKE_sgesvd_work(matrix_layout, jobu, jobvt, m, n, 
                                 storage(a), lda(a), storage(s), storage(u), lda(u),
                                 storage(vt), lda(vt), storage(work), size(work));
  }
};
  
template<>
struct gesvd_dispatch<double> {
  template<typename MATRIX, typename VECTOR>
  static lapack_int gesvd(int matrix_layout, char jobu, char jobvt,
                          lapack_int m, lapack_int n, MATRIX& a, VECTOR& s,
                          MATRIX& u, MATRIX& vt, VECTOR& work) {
    if (size(work) == std::min(m, n) - 1)
      return LAPACKE_dgesvd(matrix_layout, jobu, jobvt, m, n, 
                            storage(a), lda(a), storage(s), storage(u), lda(u),
                            storage(vt), lda(vt), storage(work));
    else
      return LAPACKE_dgesvd_work(matrix_layout, jobu, jobvt, m, n, 
                                 storage(a), lda(a), storage(s), storage(u), lda(u),
                                 storage(vt), lda(vt), storage(work), size(work));
  }
};
  
template<>
struct gesvd_dispatch<std::complex<float> > {
  template<typename MATRIX, typename VECTOR>
  static lapack_int gesvd(int matrix_layout, char jobu, char jobvt,
                          lapack_int m, lapack_int n, MATRIX& a, VECTOR& s,
                          MATRIX& u, MATRIX& vt, VECTOR& superb) {
    if (size(superb) != std::min(m, n) - 1)
      throw std::invalid_argument("vector superb size mismatch");
    return LAPACKE_cgesvd(matrix_layout, jobu, jobvt, m, n, 
                          complex_cast(storage(a)), lda(a), storage(s),
                          complex_cast(storage(u)), lda(u),
                          complex_cast(storage(vt)), lda(vt), storage(superb));
  }
  template<typename MATRIX, typename VECTOR0, typename VECTOR1>
  static lapack_int gesvd_work(int matrix_layout, char jobu, char jobvt,
                               lapack_int m, lapack_int n, MATRIX& a, VECTOR0& s,
                               MATRIX& u, MATRIX& vt, VECTOR1& work, VECTOR0& rwork) {
    return LAPACKE_cgesvd(matrix_layout, jobu, jobvt, m, n, 
                          complex_cast(storage(a)), lda(a), storage(s),
                          complex_cast(storage(u)), lda(u),
                          complex_cast(storage(vt)), lda(vt),
                          complex_cast(storage(work)), size(work), storage(rwork));
  }
};
  
template<>
struct gesvd_dispatch<std::complex<double> > {
  template<typename MATRIX, typename VECTOR>
  static lapack_int gesvd(int matrix_layout, char jobu, char jobvt,
                          lapack_int m, lapack_int n, MATRIX& a, VECTOR& s,
                          MATRIX& u, MATRIX& vt, VECTOR& superb) {
    if (size(superb) != std::min(m, n) - 1)
      throw std::invalid_argument("vector superb size mismatch");
    return LAPACKE_zgesvd(matrix_layout, jobu, jobvt, m, n, 
                          complex_cast(storage(a)), lda(a), storage(s),
                          complex_cast(storage(u)), lda(u),
                          complex_cast(storage(vt)), lda(vt), storage(superb));
  }
  template<typename MATRIX, typename VECTOR0, typename VECTOR1>
  static lapack_int gesvd_work(int matrix_layout, char jobu, char jobvt,
                               lapack_int m, lapack_int n, MATRIX& a, VECTOR0& s,
                               MATRIX& u, MATRIX& vt, VECTOR1& work, VECTOR0& rwork) {
    return LAPACKE_zgesvd(matrix_layout, jobu, jobvt, m, n, 
                          complex_cast(storage(a)), lda(a), storage(s),
                          complex_cast(storage(u)), lda(u),
                          complex_cast(storage(vt)), lda(vt),
                          complex_cast(storage(work)), size(work), storage(rwork));
  }
};
  
}
  
template<typename MATRIX, typename VECTOR>
lapack_int gesvd(char jobu, char jobvt, MATRIX& a, VECTOR& s, MATRIX& u, MATRIX& vt,
                 VECTOR& work) {
  BOOST_STATIC_ASSERT(std::is_same<typename norm_t<MATRIX>::type,
                      typename value_t<VECTOR>::type>::value);
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
  return gesvd_dispatch<typename value_t<MATRIX>::type>
    ::gesvd((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR),
            jobu, jobvt, m, n, a, s, u, vt, work);
}

template<typename MATRIX, typename VECTOR0, typename VECTOR1>
lapack_int gesvd(char jobu, char jobvt, MATRIX& a, VECTOR0& s, MATRIX& u, MATRIX& vt,
                 VECTOR1& work, VECTOR0& rwork) {
  BOOST_STATIC_ASSERT(std::is_same<typename norm_t<MATRIX>::type,
                      typename value_t<VECTOR0>::type>::value);
  BOOST_STATIC_ASSERT(std::is_same<typename value_t<MATRIX>::type,
                      typename value_t<VECTOR1>::type>::value);
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
  if (size(rwork) < 5 * r)
    throw std::invalid_argument("matrix rwork size mismatch");
  return gesvd_dispatch<typename value_t<MATRIX>::type>
    ::gesvd_work((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR),
                 jobu, jobvt, m, n, a, s, u, vt, work, rwork);
}

} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_GESVD_HPP
