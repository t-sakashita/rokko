/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_EIGEN_EXA_DIAGONALIZE_HPP
#define ROKKO_EIGEN_EXA_DIAGONALIZE_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/eigen_exa.h>

namespace rokko {
namespace eigen_exa {

template <typename MATRIX_MAJOR>
void diagonalize(rokko::distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                 rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  if(mat.is_row_major())
    throw "eigen_exa doesn't support matrix_row_major.  Use eigen_exa with matrix_col_major.";
  ROKKO_eigen_init(mat.get_grid().get_comm(), (mat.get_grid().is_row_major() ? 'R' : 'C'));
  int dim = mat.get_m_global();
  int lld = mat.get_lld();
  timer_in.start(1);
  ROKKO_eigen_sx(dim, dim, mat.get_array_pointer(), lld, &eigvals[0], eigvecs.get_array_pointer(),
                 lld, 8, 128);
  timer_in.stop(1);
  ROKKO_eigen_free(1);
}

} // namespace eigen_exa
} // namespace rokko

#endif // ROKKO_EIGEN_EXA_DIAGONALIZE_HPP
