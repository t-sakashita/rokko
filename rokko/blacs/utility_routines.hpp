/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_BLACS_UTILITY_ROUTINES_HPP
#define ROKKO_BLACS_UTILITY_ROUTINES_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/blacs/blacs_wrap.h>

namespace rokko {
namespace blacs {

template<typename MATRIX_MAJOR>
char set_grid_blacs(int ictxt, distributed_matrix<double, MATRIX_MAJOR>& mat);

template<typename MATRIX_MAJOR>
void set_desc(int ictxt, distributed_matrix<double, MATRIX_MAJOR>& mat, int desc[9]);

} // namespace blacs
} // namespace rokko

#endif // ROKKO_BLACS_UTILITY_ROUTINES_HPP

