/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_UTILITY_FRANK_MATRIX_HPP
#define ROKKO_uTILITY_FRANK_MATRIX_HPP

#include <cmath>
#include <stdexcept>
#include <boost/throw_exception.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/distributed_matrix.hpp>

namespace rokko {

class frank_matrix {
public:
  template<typename MATRIX_MAJOR>
  static void generate(rokko::localized_matrix<MATRIX_MAJOR>& mat) {
    if (mat.rows() != mat.cols())
      BOOST_THROW_EXCEPTION(std::invalid_argument("frank_matrix::generate() : non-square matrix"));
    int n = mat.rows();
    for(int i = 0; i < n; ++i) {
      for(int j = 0; j < n; ++j) {
        mat(i,j) = n - std::max(i, j);
      }
    }
  }

  template<typename MATRIX_MAJOR>
  static void generate(rokko::distributed_matrix<MATRIX_MAJOR>& mat) {
    if (mat.m_global != mat.n_global)
      BOOST_THROW_EXCEPTION(std::invalid_argument("frank_matrix::generate() : non-square matrix"));
    for(int local_i = 0; local_i < mat.m_local; ++local_i) {
      for(int local_j = 0; local_j < mat.n_local; ++local_j) {
        int global_i = mat.translate_l2g_row(local_i);
        int global_j = mat.translate_l2g_col(local_j);
        mat.set_local(local_i, local_j, mat.m_global - std::max(global_i, global_j));
      }
    }
  }

  // another (slower) implementation using set_global function
  template<typename MATRIX_MAJOR>
  static void generate_global(rokko::distributed_matrix<MATRIX_MAJOR>& mat) {
    if (mat.m_global != mat.n_global)
      BOOST_THROW_EXCEPTION(std::invalid_argument("frank_matrix::generate() : non-square matrix"));
    for(int global_i=0; global_i<mat.m_global; ++global_i) {
      for(int global_j=0; global_j<mat.n_global; ++global_j) {
        mat.set_global(global_i, global_j, mat.m_global - std::max(global_i, global_j) );
      }
    }
  }

  // calculate k-th smallest eigenvalue of dim-dimensional Frank matrix (k=0...dim-1)
  static double eigenvalue(int dim, int k) {
    return 1 / (2 * (1 - std::cos(M_PI * (2 * k + 1) / (2 * dim + 1))));
  }
};
    
} // namespace rokko

#endif // ROKKO_UTILITY_FRANK_MATRIX_HPP
