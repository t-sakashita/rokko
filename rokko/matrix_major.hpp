/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MATRIX_MAJOR_HPP
#define ROKKO_MATRIX_MAJOR_HPP

namespace rokko {

extern struct matrix_row_major {} matrix_row_major_d;
extern struct matrix_col_major {} matrix_col_major_d;

} // namespace rokko

#endif // ROKKO_MATRIX_MAJOR_HPP
