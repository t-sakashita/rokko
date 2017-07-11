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

#ifndef ROKKO_LAPACK_COMPLEX_CAST_HPP
#define ROKKO_LAPACK_COMPLEX_CAST_HPP

#include <complex>

namespace rokko {
namespace lapack {

#define COMPLEX_CAST_IMPL(T) \
inline T * complex_cast(T* data) { return data; } \
inline const T * complex_cast(const T* data) { return data; } \
inline lapack_complex_ ## T * complex_cast(std::complex< T >* data) { return reinterpret_cast<lapack_complex_ ## T *>(data); } \
inline const lapack_complex_ ## T * complex_cast(const std::complex< T >* data) { return reinterpret_cast<const lapack_complex_ ## T *>(data); } \
inline const T & complex_cast(const T & data) { return *complex_cast(&data); } \
inline const lapack_complex_ ## T & complex_cast(const std::complex< T >& data) { return *complex_cast(&data); }

COMPLEX_CAST_IMPL(float)
COMPLEX_CAST_IMPL(double)

} // end namespace lapack
} // end namespace rokko

#endif // ROKKO_LAPACK_COMPLEX_CAST_HPP
