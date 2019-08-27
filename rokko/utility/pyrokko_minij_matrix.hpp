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

#ifndef PYROKKO_MINIJ_MATRIX_HPP
#define PYROKKO_MINIJ_MATRIX_HPP

#include <rokko/utility/minij_matrix.hpp>

#include <rokko/pyrokko_localized_matrix.hpp>
#include <rokko/pyrokko_distributed_matrix.hpp>

namespace rokko {

class wrap_minij_matrix : public minij_matrix {
public:
  static void generate(wrap_localized_matrix& mat) {
    if (mat.is_major_col())
      minij_matrix::generate(mat.col_ver());
    else
      minij_matrix::generate(mat.row_ver());
  }

  static void generate(wrap_distributed_matrix& mat) {
    if (mat.is_major_col())
      minij_matrix::generate(mat.col_ver());
    else
      minij_matrix::generate(mat.row_ver());
  }};
} // namespace rokko

#endif // PYROKKO_MINIJ_MATRIX_HPP
