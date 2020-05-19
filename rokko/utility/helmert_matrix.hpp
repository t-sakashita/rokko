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

#ifndef ROKKO_UTILITY_HELMERT_MATRIX_HPP
#define ROKKO_UTILITY_HELMERT_MATRIX_HPP

#include <cmath>
#include <stdexcept>
#include <rokko/config.h>
#include <rokko/eigen3.hpp>
#if defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
# include <rokko/distributed_matrix.hpp>
#endif

namespace rokko {

class helmert_matrix {
public:
  template<typename T, int ROWS, int COLS, int MATRIX_MAJOR>
  static void generate(Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& mat) {
    if (mat.rows() != mat.cols())
      throw std::invalid_argument("helmert_matrix::generate() : non-square matrix");
    const int n = mat.rows();
    mat.setZero();
    mat.row(0).fill( 1 / sqrt(n) );
    for (int i=1; i < mat.rows(); ++i) {
      mat.row(i).head(i).fill( 1 / sqrt(static_cast<T>(i*(i+1))) );
      mat(i,i) = - sqrt(static_cast<T>(i)/(i+1));
    }
  }

  template<typename T, int ROWS, int COLS, int MATRIX_MAJOR, int SIZE>
  static void generate_for_given_eigenvalues(Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& mat, Eigen::Vector<T, SIZE> const& diag) {
    if (mat.rows() != mat.cols())
      throw std::invalid_argument("helmert_matrix::generate_for_given_eigenvalues() : non-square matrix");
    const int n = mat.rows();
    for (int i=0; i<n; ++i) {
      T common_elem = diag(0) / n;
      for (int k=i+1; k<n; ++k)
        common_elem += diag(k) / (k*(k+1));  // Remark: i=max(i,j)
      mat(i, i) = common_elem + i * diag(i) / (i+1);
      T val = common_elem - diag(i) / (i+1);
      mat.row(i).head(i).setConstant(val);
      mat.col(i).head(i).setConstant(val);
    }
  }
  
#if defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
  template<typename T, typename MATRIX_MAJOR>
  static void generate(rokko::distributed_matrix<T, MATRIX_MAJOR>& mat) {
    if (mat.get_m_global() != mat.get_n_global())
      throw std::invalid_argument("helmert_matrix::generate() : non-square matrix");
    const int n = mat.get_m_global();
    int start_i;
    if (mat.has_global_row_index(0)) {
      T val = 1 / sqrt(n);
      for(int local_j = 0; local_j < mat.get_n_local(); ++local_j)
        mat.set_local(0, local_j, val);
      start_i = 1;
    }
    else start_i = 0;

    for(int local_i = start_i; local_i < mat.get_m_local(); ++local_i) {
      int global_i = mat.translate_l2g_row(local_i);
      T val = 1 / sqrt(static_cast<T>(global_i*(global_i+1)));
      for(int local_j = 0; local_j < mat.get_n_local(); ++local_j) {
        int global_j = mat.translate_l2g_col(local_j);
	if (global_j < global_i) mat.set_local(local_i, local_j, val);
	else if (global_j == global_i) mat.set_local(local_i, local_j, - sqrt(static_cast<T>(global_i)/(global_i+1)));
      }
    }
  }

  template<typename T, typename MATRIX_MAJOR, int SIZE>
  static void generate_for_given_eigenvalues(rokko::distributed_matrix<T, MATRIX_MAJOR>& mat, Eigen::Vector<T, SIZE> const& diag) {
    if (mat.get_m_global() != mat.get_n_global())
      throw std::invalid_argument("helmert_matrix::generate_for_given_eigenvalues() : non-square matrix");
    const int n = mat.get_m_global();

    for(int local_i = 0; local_i < mat.get_m_local(); ++local_i) {
      int global_i = mat.translate_l2g_row(local_i);
      T common_elem = diag(0) / n;
      for (int k=global_i+1; k<n; ++k)
	common_elem += diag(k) / (k*(k+1));  // Remark: i=max(i,j)
      T val = common_elem - diag(global_i) / (global_i+1);
      for(int local_j = 0; local_j < mat.get_n_local(); ++local_j) {
        int global_j = mat.translate_l2g_col(local_j);
	if (global_j < global_i) mat.set_local(local_i, local_j, val);
	else if (global_j == global_i) mat.set_local(local_i, local_j, common_elem + global_i * diag(global_i) / (global_i+1));
      }
    }
    
    for(int local_j = 0; local_j < mat.get_n_local(); ++local_j) {
      int global_j = mat.translate_l2g_col(local_j);
      T common_elem = diag(0) / n;
      for (int k=global_j+1; k<n; ++k)
	common_elem += diag(k) / (k*(k+1));  // Remark: i=max(i,j)
      T val = common_elem - diag(global_j) / (global_j+1);
      for(int local_i = 0; local_i < mat.get_m_local(); ++local_i) {
        int global_i = mat.translate_l2g_row(local_i);
	if (global_i < global_j) mat.set_local(local_i, local_j, val);
      }
    }
  }
#endif
};
    
} // namespace rokko

#endif // ROKKO_UTILITY_HELMERT_MATRIX_HPP
