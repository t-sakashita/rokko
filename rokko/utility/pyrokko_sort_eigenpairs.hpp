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

#ifndef PYROKKO_UTILITY_SORT_EIGENPAIRS_HPP
#define PYROKKO_UTILITY_SORT_EIGENPAIRS_HPP

#include <rokko/utility/sort_eigenpairs.hpp>

#include <rokko/pyrokko_localized_matrix.hpp>

namespace rokko {

void pyrokko_sort_eigenpairs(const wrap_localized_vector& eigval,
                             const wrap_localized_matrix& eigvec,
                             wrap_localized_vector& eigval_sorted,
                             wrap_localized_matrix& eigvec_sorted,
                             bool ascending = true) {
  if (eigvec.is_major_col())
    sort_eigenpairs(eigval.obj(), eigvec.col_ver(), eigval_sorted.obj(), eigvec_sorted.col_ver(), ascending);
  else
    sort_eigenpairs(eigval.obj(), eigvec.row_ver(), eigval_sorted.obj(), eigvec_sorted.row_ver(), ascending);
}

} // namespace rokko

#endif // PYROKKO_UTILITY_SORT_EIGENPAIRS_HPP
