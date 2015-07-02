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

#ifndef ROKKO_UTILITY_HELMERT_MATRIX_HPP
#define ROKKO_UTILITY_HELMERT_MATRIX_HPP

#include <cmath>
#include <stdexcept>
#include <boost/throw_exception.hpp>
#include <rokko/config.h>
#include <rokko/localized_matrix.hpp>
#if defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
# include <rokko/distributed_matrix.hpp>
#endif

namespace rokko {

class helmert_matrix {
public:
  template<typename T, typename MATRIX_MAJOR>
  static void generate(rokko::localized_matrix<T, MATRIX_MAJOR>& mat) {
    if (mat.rows() != mat.cols())
      BOOST_THROW_EXCEPTION(std::invalid_argument("helmert_matrix::generate() : non-square matrix"));
    int n = mat.rows();
    mat.row(0).fill( 1 / sqrt(n) );
    for (int i=1; i < mat.rows(); ++i) {
      mat.row(i).head(i).fill( 1 / sqrt(static_cast<T>(i*(i+1))) );
      mat(i,i) = - sqrt(static_cast<T>(i)/(i+1));
    }
  }

  template<typename T, typename MATRIX_MAJOR>
  static void generate_for_given_eigenvalues(rokko::localized_matrix<T, MATRIX_MAJOR>& mat, rokko::localized_vector<T> const& diag) {
    if (mat.rows() != mat.cols())
      BOOST_THROW_EXCEPTION(std::invalid_argument("helmert_matrix::generate() : non-square matrix"));
    int n = mat.rows();
    for (int i=0; i<n; ++i) {
      double common_elem = diag(0) / n;
      for (int k=i+1; k<n; ++k)
	common_elem += diag(k) / (k*(k+1));  // Remark: i=max(i,j)
      mat(i, i) = common_elem + i * diag(i) / (i+1);
      double val = common_elem - diag(i) / (i+1);
      mat.row(i).head(i).setConstant(val);
      mat.col(i).head(i).setConstant(val);
    }
  }
  
#if defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
  template<typename T, typename MATRIX_MAJOR>
  static void generate(rokko::distributed_matrix<T, MATRIX_MAJOR>& mat) {
    if (mat.get_m_global() != mat.get_n_global())
      BOOST_THROW_EXCEPTION(std::invalid_argument("helmert_matrix::generate() : non-square matrix"));
    const int n = mat.get_m_global();
    int start_i;
    if (mat.is_gindex_myrow(0)) {
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

  /*
  // another (slower) implementation using set_global function
  template<typename T, typename MATRIX_MAJOR>
  static void generate_global(rokko::distributed_matrix<MATRIX_MAJOR>& mat) {
    if (mat.m_global != mat.n_global)
      BOOST_THROW_EXCEPTION(std::invalid_argument("helmert_matrix::generate() : non-square matrix"));
    for(int global_i=0; global_i<mat.m_global; ++global_i) {
      for(int global_j=0; global_j<mat.n_global; ++global_j) {
        mat.set_global(global_i, global_j, mat.m_global - std::max(global_i, global_j) );
      }
    }
  }
  */

  template<typename T, typename MATRIX_MAJOR>
  static void generate_for_given_eigenvalues(rokko::distributed_matrix<T, MATRIX_MAJOR>& mat, rokko::localized_vector<T> const& diag) {
    if (mat.get_m_global() != mat.get_n_global())
      BOOST_THROW_EXCEPTION(std::invalid_argument("helmert_matrix::generate() : non-square matrix"));
    const int n = mat.get_m_global();

    for(int local_i = 0; local_i < mat.get_m_local(); ++local_i) {
      int global_i = mat.translate_l2g_row(local_i);
      double common_elem = diag(0) / n;
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
      double common_elem = diag(0) / n;
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
