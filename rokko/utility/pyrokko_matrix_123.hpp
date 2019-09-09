/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef PYROKKO_MATRIX_123_HPP
#define PYROKKO_MATRIX_123_HPP

#include <rokko/utility/matrix123.hpp>

#include <rokko/pyrokko_distributed_matrix.hpp>

namespace rokko {

static void pyrokko_generate_matrix_123(wrap_distributed_matrix& mat) {
  if (mat.is_major_col())
    generate_matrix_123(mat.col_ver());
  else
    generate_matrix_123(mat.row_ver());
}

} // namespace rokko

#endif // PYROKKO_MATRIX_123_HPP
