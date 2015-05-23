/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2015 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_LAPACK_CONFIG_H
#define ROKKO_LAPACK_CONFIG_H

#ifdef __cplusplus
# include <complex>
# ifndef lapack_complex_float
#   define lapack_complex_float std::complex<float>
# endif
# ifndef lapack_complex_double
#   define lapack_complex_double std::complex<double>
# endif
#endif

#include <lapacke_mangling.h>

#endif // ROKKO_LAPACK_CONFIG_H
