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

#ifndef ROKKO_UTILITY_MATRIX123_HPP
#define ROKKO_UTILITY_MATRIX123_HPP

namespace rokko {

template<typename T>
void generate_matrix_123(rokko::distributed_matrix<T>& mat) {
  for (int local_i = 0; local_i < mat.m_local; ++local_i) {
    for (int local_j = 0; local_j < mat.n_local; ++local_j) {
      int global_i = mat.translate_l2g_row(local_i);
      int global_j = mat.translate_l2g_col(local_j);
      mat.set_local(local_i, local_j, mat.m_global * global_j + global_i );
    }
  }
}

} // namespace rokko

#endif // ROKKO_UTILITY_MATRIX123_HPP
