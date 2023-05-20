/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <cmath>
#include <stdexcept>
#include <rokko/config.h>
#include <rokko/eigen3.hpp>
#if defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
# include <rokko/distributed_matrix.hpp>
#endif

namespace rokko {

class minij_matrix {
public:
  template<typename T, int ROWS, int COLS>
  static void generate(Eigen::Matrix<T,ROWS,COLS,Eigen::ColMajor>& mat) {
    if (mat.rows() != mat.cols())
      throw std::invalid_argument("minij_matrix::generate() : non-square matrix");
    const auto n = mat.rows();
    for(int j = 0; j < n; ++j) {
      for(int i = 0; i < n; ++i) {
        mat(i,j) = std::min(i, j) + 1;
      }
    }
  }

  template<typename T, int ROWS, int COLS>
  static void generate(Eigen::Matrix<T,ROWS,COLS,Eigen::RowMajor>& mat) {
    if (mat.rows() != mat.cols())
      throw std::invalid_argument("minij_matrix::generate() : non-square matrix");
    const auto n = mat.rows();
    for(int i = 0; i < n; ++i) {
      for(int j = 0; j < n; ++j) {
        mat(i,j) = std::min(i, j) + 1;
      }
    }
  }

#if defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
  template<typename T>
  static void generate(rokko::distributed_matrix<T, matrix_col_major>& mat) {
    auto const& map = mat.get_mapping();
    if (map.get_m_global() != map.get_n_global())
      throw std::invalid_argument("minij_matrix::generate() : non-square matrix");

    for(int local_j = 0; local_j < map.get_n_local(); ++local_j) {
      int global_j = map.translate_l2g_col(local_j);
      for(int local_i = 0; local_i < map.get_m_local(); ++local_i) {
        int global_i = map.translate_l2g_row(local_i);
        mat.set_local(local_i, local_j, std::min(global_i, global_j) + 1);
      }
    }
  }

  template<typename T>
  static void generate(rokko::distributed_matrix<T, matrix_row_major>& mat) {
    auto const& map = mat.get_mapping();
    if (map.get_m_global() != map.get_n_global())
      throw std::invalid_argument("minij_matrix::generate() : non-square matrix");

    for(int local_i = 0; local_i < map.get_m_local(); ++local_i) {
      int global_i = map.translate_l2g_row(local_i);
      for(int local_j = 0; local_j < map.get_n_local(); ++local_j) {
        int global_j = map.translate_l2g_col(local_j);
        mat.set_local(local_i, local_j, std::min(global_i, global_j) + 1);
      }
    }
  }

  // another (slower) implementation using set_global function
  template<typename T>
  static void generate_global(rokko::distributed_matrix<T, matrix_col_major>& mat) {
    auto const& map = mat.get_mapping();
    if (map.m_global() != map.n_global())
      throw std::invalid_argument("minij_matrix::generate() : non-square matrix");

    for(int global_j=0; global_j<map.n_global(); ++global_j) {
      for(int global_i=0; global_i<map.m_global(); ++global_i) {
        mat.set_global(global_i, global_j, std::min(global_i, global_j) + 1);
      }
    }
  }

  template<typename T>
  static void generate_global(rokko::distributed_matrix<T, matrix_row_major>& mat) {
    auto const& map = mat.get_mapping();
    if (map.m_global() != map.n_global())
      throw std::invalid_argument("minij_matrix::generate() : non-square matrix");

    for(int global_i=0; global_i<map.m_global(); ++global_i) {
      for(int global_j=0; global_j<map.n_global(); ++global_j) {
        mat.set_global(global_i, global_j, std::min(global_i, global_j) + 1);
      }
    }
  }
#endif
  
  // calculate k-th smallest eigenvalue of dim-dimensional Minij matrix (k=0...dim-1)
  static double eigenvalue(int dim, int k) {
    return 1 / (2 * (1 - std::cos(M_PI * (2 * k + 1) / (2 * dim + 1))));
  }
};

} // namespace rokko
