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

#ifndef ROKKO_EIGENEXA_DIAGONALIZE_EIGEN_S_HPP
#define ROKKO_EIGENEXA_DIAGONALIZE_EIGEN_S_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/eigenexa.hpp>
#include <rokko/utility/timer.hpp>

namespace rokko {
namespace eigenexa {

// eigen_s eigenvalues / eigenvectors
template <typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_eigen_s(rokko::distributed_matrix<double, MATRIX_MAJOR>& mat,
			       VEC& eigvals, rokko::distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			       parameters const& params) {
  parameters params_out;
  if(mat.is_row_major())
    throw std::invalid_argument("eigenexa::diagonalize_eigen_s() : eigenexa doesn't support matrix_row_major.  Use eigenexa with matrix_col_major.");
  if((mat.get_mb() != 1) || (mat.get_nb() != 1))
    throw std::invalid_argument("eigenexa::diagonalize_eigen_s() : eigenexa supports only 1x1 block size.");
  rokko::eigenexa::init(mat.get_grid().get_comm(), (mat.get_grid().is_row_major() ? 'R' : 'C'));
  int m_forward = 48, m_backward = 128;
  get_key(params, "m_forward", m_forward);
  get_key(params, "m_backward", m_backward);

  rokko::eigenexa::eigen_s(mat, eigvals, eigvecs, m_forward, m_backward, 'A');

  rokko::eigenexa::free(1);
  return params_out;
}

// eigen_s only eigenvalues
template <typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_eigen_s(rokko::distributed_matrix<double, MATRIX_MAJOR>& mat,
			       VEC& eigvals,
			       parameters const& params) {
  parameters params_out;
  if(mat.is_row_major())
    throw std::invalid_argument("eigenexa::diagonalize_eigen_s() : eigenexa doesn't support matrix_row_major.  Use eigenexa with matrix_col_major.");
  if((mat.get_mb() != 1) || (mat.get_nb() != 1))
    throw std::invalid_argument("eigenexa::diagonalize_eigen_s() : eigenexa supports only 1x1 block size.");
  rokko::eigenexa::init(mat.get_grid().get_comm(), (mat.get_grid().is_row_major() ? 'R' : 'C'));
  int m_forward = 48, m_backward = 128;
  get_key(params, "m_forward", m_forward);
  get_key(params, "m_backward", m_backward);

  rokko::eigenexa::eigen_s(mat, eigvals, m_forward, m_backward, 'N');

  rokko::eigenexa::free(1);
  return params_out;
}

} // namespace eigenexa
} // namespace rokko

#endif // ROKKO_EIGENEXA_DIAGONALIZE_EIGEN_S_HPP
