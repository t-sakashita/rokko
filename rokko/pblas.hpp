/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2017 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_PBLAS_HPP
#define ROKKO_PBLAS_HPP

#include <rokko/cpblas.h>
#include <rokko/lapack/complex_cast.hpp>

using rokko::lapack::complex_cast;

namespace rokko {
namespace pblas {

// pcopy

namespace {

#define PBLAS_PCOPY_IMPL(NAMES, TYPE)                                   \
inline void pcopy_dispatch(int n, const TYPE * x, int ix, int jx, const int* descx, int incx, TYPE * y, int iy, int jy, const int* descy, int incy) { cpblas_ ## NAMES (n, complex_cast(x), ix, jx, descx, incx, complex_cast(y), iy, jy, descy, incy); }
    
PBLAS_PCOPY_IMPL(pscopy, float);
PBLAS_PCOPY_IMPL(pdcopy, double);
PBLAS_PCOPY_IMPL(pccopy, std::complex<float>);
PBLAS_PCOPY_IMPL(pzcopy, std::complex<double>);

#undef PBLAS_PCOPY_IMPL

}

template<typename MATRIX>
void pcopy(int n, const MATRIX& x, int ix, int jx, int incx, MATRIX& y, int iy, int jy, int incy) {
  const int* descx = x.get_mapping().get_blacs_descriptor();
  const int* descy = y.get_mapping().get_blacs_descriptor();
  pcopy_dispatch(n, x.get_array_pointer(), ix, jx, descx, incx,
                 y.get_array_pointer(), iy, jy, descy, incy);
}

// pdot, pdotu, pdotc

namespace {
    
#define PBLAS_PDOT_IMPL(NAMEX, NAMES, TYPE) \
inline TYPE NAMEX ## _dispatch (int n, const TYPE * x, int ix, int jx, const int* descx, int incx, const TYPE * y, int iy, int jy, const int* descy, int incy) { \
  TYPE dot; \
  cpblas_ ## NAMES ## _sub (n, complex_cast(&dot), complex_cast(x), ix, jx, descx, incx, complex_cast(y), iy, jy, descy, incy); \
  return dot; \
}

PBLAS_PDOT_IMPL(pdot, psdot, float);
PBLAS_PDOT_IMPL(pdot, pddot, double);
PBLAS_PDOT_IMPL(pdot, pcdotc, std::complex<float>);
PBLAS_PDOT_IMPL(pdot, pzdotc, std::complex<double>);
PBLAS_PDOT_IMPL(pdotu, pcdotu, std::complex<float>);
PBLAS_PDOT_IMPL(pdotu, pzdotu, std::complex<double>);
PBLAS_PDOT_IMPL(pdotc, psdot, float);
PBLAS_PDOT_IMPL(pdotc, pddot, double);
PBLAS_PDOT_IMPL(pdotc, pcdotc, std::complex<float>);
PBLAS_PDOT_IMPL(pdotc, pzdotc, std::complex<double>);

#undef PBLAS_PDOT_IMPL

}
  
template<typename VECTOR>
typename VECTOR::value_type pdot(int n, const VECTOR& x, int ix, int jx, int incx,
                                 const VECTOR& y, int iy, int jy, int incy) {
  const int* descx = x.get_mapping().get_blacs_descriptor();
  const int* descy = y.get_mapping().get_blacs_descriptor();
  return pdot_dispatch(n, x.get_array_pointer(), ix, jx, descx, incx,
                       y.get_array_pointer(), iy, jy, descy, incy);
}

template<typename VECTOR>
typename VECTOR::value_type pdotc(int n, const VECTOR& x, int ix, int jx, int incx,
                                 const VECTOR& y, int iy, int jy, int incy) {
  const int* descx = x.get_mapping().get_blacs_descriptor();
  const int* descy = y.get_mapping().get_blacs_descriptor();
  return pdotc_dispatch(n, x.get_array_pointer(), ix, jx, descx, incx,
                        y.get_array_pointer(), iy, jy, descy, incy);
}

template<typename VECTOR>
typename VECTOR::value_type pdotu(int n, const VECTOR& x, int ix, int jx, int incx,
                                 const VECTOR& y, int iy, int jy, int incy) {
  const int* descx = x.get_mapping().get_blacs_descriptor();
  const int* descy = y.get_mapping().get_blacs_descriptor();
  return pdotu_dispatch(n, x.get_array_pointer(), ix, jx, descx, incx,
                        y.get_array_pointer(), iy, jy, descy, incy);
}

// pgemv

namespace {
    
#define PBLAS_PGEMV_IMPL(NAMES, TYPE) \
inline void pgemv_dispatch(char trans, int m, int n, TYPE alpha, const TYPE* a, int ia, int ja, const int* desca, const TYPE* x, int ix, int jx, const int* descx, int incx, TYPE beta, TYPE * y, int iy, int jy, const int* descy, int incy) { \
  cpblas_ ## NAMES (trans, m, n, complex_cast(alpha), complex_cast(a), ia, ja, desca, complex_cast(x), ix, jx, descx, incx, complex_cast(beta), complex_cast(y), iy, jy, descy, incy); \
}

PBLAS_PGEMV_IMPL(psgemv, float);
PBLAS_PGEMV_IMPL(pdgemv, double);
PBLAS_PGEMV_IMPL(pcgemv, std::complex<float>);
PBLAS_PGEMV_IMPL(pzgemv, std::complex<double>);

#undef PBLAS_PGEMV_IMPL

}
  
template<typename MATRIX, typename VECTOR, typename T>
void pgemv(char trans, T alpha, const MATRIX& a, const VECTOR& x, int incx,
           T beta, VECTOR& y, int incy) {
  const int* desca = a.get_mapping().get_blacs_descriptor();
  const int* descx = x.get_mapping().get_blacs_descriptor();
  const int* descy = y.get_mapping().get_blacs_descriptor();
  pgemv_dispatch(trans, a.get_m_global(), a.get_n_global(), alpha,
                 a.get_array_pointer(), 0, 0, desca,
                 x.get_array_pointer(), 0, 0, descx, incx, beta,
                 y.get_array_pointer(), 0, 0, descy, incy);
}

// pgemm

namespace {
    
#define PBLAS_PGEMM_IMPL(NAMES, TYPE) \
inline void pgemm_dispatch(char transa, char transb, int m, int n, int k, TYPE alpha, const TYPE* a, int ia, int ja, const int* desca, const TYPE* b, int ib, int jb, const int* descb, TYPE beta, TYPE* c, int ic, int jc, const int* descc) { \
  cpblas_ ## NAMES (transa, transb, m, n, k, complex_cast(alpha), complex_cast(a), ia, ja, desca, complex_cast(b), ib, jb, descb, complex_cast(beta), complex_cast(c), ic, jc, descc); \
}

PBLAS_PGEMM_IMPL(psgemm, float);
PBLAS_PGEMM_IMPL(pdgemm, double);
PBLAS_PGEMM_IMPL(pcgemm, std::complex<float>);
PBLAS_PGEMM_IMPL(pzgemm, std::complex<double>);

#undef PBLAS_PGEMM_IMPL

}
  
template<typename MATRIX, typename T>
void pgemm(char transa, char transb, T alpha, const MATRIX& a, const MATRIX& b,
           T beta, MATRIX& c) {
  const int* desca = a.get_mapping().get_blacs_descriptor();
  const int* descb = b.get_mapping().get_blacs_descriptor();
  const int* descc = c.get_mapping().get_blacs_descriptor();
  pgemm_dispatch(transa, transb, a.get_m_global(), b.get_n_global(), a.get_n_global(), alpha,
                 a.get_array_pointer(), 0, 0, desca,
                 b.get_array_pointer(), 0, 0, descb, beta,
                 c.get_array_pointer(), 0, 0, descc);
}

} // namespace pblas
} // namespace rokko

#endif // ROKKO_PBLAS_HPP
