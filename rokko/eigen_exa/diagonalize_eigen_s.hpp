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

#ifndef ROKKO_EIGEN_EXA_DIAGONALIZE_EIGEN_S_HPP
#define ROKKO_EIGEN_EXA_DIAGONALIZE_EIGEN_S_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/eigen_exa/eigen_exa_wrap.h>
#include <rokko/utility/timer.hpp>

namespace rokko {
namespace eigen_exa {

// eigen_s eigenvalues / eigenvectors
template <typename MATRIX_MAJOR>
void diagonalize_eigen_s(rokko::distributed_matrix<double, MATRIX_MAJOR>& mat,
			 localized_vector<double>& eigvals, rokko::distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			 parameters const& params, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  if(mat.is_row_major())
    throw "eigen_exa doesn't support matrix_row_major.  Use eigen_exa with matrix_col_major.";
  ROKKO_eigen_exa_init(mat.get_grid().get_comm(), (mat.get_grid().is_row_major() ? 'R' : 'C'));
  int dim = mat.get_m_global();
  int lld = mat.get_lld();
  int m_forward = 48, m_backward = 128;
  get_key(params, "m_forward", m_forward);
  get_key(params, "m_backward", m_backward);
  timer.stop(timer_id::diagonalize_initialize);
  timer.start(timer_id::diagonalize_diagonalize);
  ROKKO_eigen_exa_s(dim, dim, mat.get_array_pointer(), lld, &eigvals[0],
                    eigvecs.get_array_pointer(), lld, m_forward, m_backward, 'A');
  timer.stop(timer_id::diagonalize_diagonalize);
  timer.start(timer_id::diagonalize_finalize);
  ROKKO_eigen_exa_free(1);
  timer.stop(timer_id::diagonalize_finalize);
}

// eigen_s only eigenvalues
template <typename MATRIX_MAJOR>
void diagonalize_eigen_s(rokko::distributed_matrix<double, MATRIX_MAJOR>& mat,
			 localized_vector<double>& eigvals,
			 parameters const& params, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  if(mat.is_row_major())
    throw "eigen_exa doesn't support matrix_row_major.  Use eigen_exa with matrix_col_major.";
  ROKKO_eigen_exa_init(mat.get_grid().get_comm(), (mat.get_grid().is_row_major() ? 'R' : 'C'));
  int dim = mat.get_m_global();
  int lld = mat.get_lld();
  int m_forward = 48, m_backward = 128;
  get_key(params, "m_forward", m_forward);
  get_key(params, "m_backward", m_backward);
  timer.stop(timer_id::diagonalize_initialize);
  timer.start(timer_id::diagonalize_diagonalize);
  ROKKO_eigen_exa_s(dim, dim, mat.get_array_pointer(), lld, &eigvals[0], NULL,
		    lld, m_forward, m_backward, 'N');
  timer.stop(timer_id::diagonalize_diagonalize);
  timer.start(timer_id::diagonalize_finalize);
  ROKKO_eigen_exa_free(1);
  timer.stop(timer_id::diagonalize_finalize);
}

} // namespace eigen_exa
} // namespace rokko

#endif // ROKKO_EIGEN_EXA_DIAGONALIZE_EIGEN_S_HPP
