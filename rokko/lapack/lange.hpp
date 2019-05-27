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

#ifndef ROKKO_LAPACK_LANGE_HPP
#define ROKKO_LAPACK_LANGE_HPP

#include <complex>
#include <stdexcept>
#include <lapacke.h>
#undef I
#include <rokko/traits/norm_t.hpp>
#include "complex_cast.hpp"

namespace rokko {
namespace lapack {

namespace {

template<typename T> struct lange_dispatch;
  
template<>
struct lange_dispatch<float> {
  template<typename MATRIX>
  static float lange(int matrix_layout, char norm, lapack_int m, lapack_int n,
                     MATRIX const& a) {
    return LAPACKE_slange(matrix_layout, norm, m, n, storage(a), lda(a));
  }
};

template<>
struct lange_dispatch<double> {
  template<typename MATRIX>
  static double lange(int matrix_layout, char norm, lapack_int m, lapack_int n,
                      MATRIX const& a) {
    return LAPACKE_dlange(matrix_layout, norm, m, n, storage(a), lda(a));
  }
};

template<>
struct lange_dispatch<std::complex<float> > {
  template<typename MATRIX>
  static float lange(int matrix_layout, char norm, lapack_int m, lapack_int n,
                     MATRIX const& a) {
    return LAPACKE_clange(matrix_layout, norm, m, n, complex_cast((storage(a))), lda(a));
  }
};

template<>
struct lange_dispatch<std::complex<double> > {
  template<typename MATRIX>
  static double lange(int matrix_layout, char norm, lapack_int m, lapack_int n,
                      MATRIX const& a) {
    return LAPACKE_zlange(matrix_layout, norm, m, n, complex_cast((storage(a))), lda(a));
  }
};

}

template<typename MATRIX>
typename norm_t<MATRIX>::type lange(char norm, MATRIX const& a) {
  lapack_int m = rows(a);
  lapack_int n = cols(a);
  return lange_dispatch<typename value_t<MATRIX>::type>
    ::lange((is_col_major(a) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR), norm,
            m, n, a);
}

} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_LANGE_HPP
