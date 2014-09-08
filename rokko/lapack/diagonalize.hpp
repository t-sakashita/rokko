/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_LAPACK_DIAGONALIZE_H
#define ROKKO_LAPACK_DIAGONALIZE_H

#include "rokko/localized_matrix.hpp"
#include "rokko/localized_vector.hpp"
#include <lapacke.h>

namespace rokko {
namespace lapack {

template<typename MATRIX_MAJOR>
int diagonalize(localized_matrix<MATRIX_MAJOR>& mat, double* eigvals,
                localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  int info;
  int dim = mat.rows();
  timer_in.start(1);
  if(mat.is_col_major())
    info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', dim, &mat(0,0), dim, eigvals);
  else
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', dim, &mat(0,0), dim, eigvals);
  timer_in.stop(1);
  eigvecs = mat;
  if (info) {
    std::cerr << "error at dsyev function. info=" << info  << std::endl;
    exit(1);
  }
  return info;
}

template<class MATRIX_MAJOR>
int diagonalize(localized_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  int dim = mat.rows();
  if (eigvals.size() < dim) eigvals.resize(dim);
  return diagonalize(mat, &eigvals[0], eigvecs, timer_in);
}

template<class MATRIX_MAJOR>
int diagonalize(localized_matrix<MATRIX_MAJOR>& mat, std::vector<double>& eigvals,
                localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  int dim = mat.rows();
  if (eigvals.size() < dim) eigvals.resize(dim);
  return diagonalize(mat, &eigvals[0], eigvecs, timer_in);
}

} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_H
