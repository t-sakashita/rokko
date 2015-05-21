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

#ifndef ROKKO_LAPACKE_HPP
#define ROKKO_LAPACKE_HPP

#include <complex>
#ifndef lapack_complex_float
# define lapack_complex_float std::complex<float>
#endif
#ifndef lapack_complex_double
# define lapack_complex_double std::complex<double>
#endif

#include <lapacke.h>

#endif // ROKKO_LAPACKE_HPP
