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

#ifndef PYROKKO_COLLECTIVE_HPP
#define PYROKKO_COLLECTIVE_HPP

#include <rokko/pyrokko_distributed_matrix.hpp>
#include <rokko/pyrokko_localized_matrix.hpp>

#include <rokko/collective.hpp>

namespace rokko {

void pyrokko_gather(wrap_distributed_matrix const& from, wrap_localized_matrix& to, int root) {
  assert(from.is_major_col() == to.is_major_col());
  bool is_col = from.is_major_col();
  
  if (is_col)
    gather(from.col_ver(), to.col_ver(), root);
  else
    gather(from.row_ver(), to.row_ver(), root);
}

void pyrokko_scatter(wrap_localized_matrix const& from, wrap_distributed_matrix& to, int root) {
  assert(from.is_major_col() == to.is_major_col());
  bool is_col = from.is_major_col();
  
  if (is_col)
    scatter(from.col_ver(), to.col_ver(), root);
  else
    scatter(from.row_ver(), to.row_ver(), root);
}

} // end namespace rokko

#endif // PYROKKO_COLLECTIVE_HPP
