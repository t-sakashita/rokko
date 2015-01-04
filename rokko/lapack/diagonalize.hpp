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

#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>
#include <lapacke.h>

namespace rokko {
namespace lapack {

template<typename MATRIX_MAJOR>
int diagonalize(localized_matrix<MATRIX_MAJOR>& mat, double* eigvals,
                localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer) {
  timer.start(timer_id::diagonalize_diagonalize);
  int dim = mat.rows();
  int info;
  if(mat.is_col_major())
    info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', dim, &mat(0,0), dim, eigvals);
  else
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', dim, &mat(0,0), dim, eigvals);
  timer.stop(timer_id::diagonalize_diagonalize);
  timer.start(timer_id::diagonalize_finalize);
  eigvecs = mat;
  if (info) {
    std::cerr << "error at dsyev function. info=" << info  << std::endl;
    exit(1);
  }
  timer.stop(timer_id::diagonalize_finalize);
  return info;
}

template<typename MATRIX_MAJOR, typename VEC>
int diagonalize(localized_matrix<MATRIX_MAJOR>& mat, VEC& eigvals,
                localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  int dim = mat.rows();
  if (eigvals.size() < dim) eigvals.resize(dim);
  timer.stop(timer_id::diagonalize_initialize);
  return diagonalize(mat, &eigvals[0], eigvecs, timer);
}

} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_H
