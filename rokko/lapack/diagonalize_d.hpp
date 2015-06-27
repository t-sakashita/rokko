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

#ifndef ROKKO_LAPACK_DIAGONALIZE_D_HPP
#define ROKKO_LAPACK_DIAGONALIZE_D_HPP

#include <rokko/parameters.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.h>

namespace rokko {
namespace lapack {

// dsyevd only eigenvalues
template<typename MATRIX_MAJOR>
int diagonalize_d(localized_matrix<double, MATRIX_MAJOR>& mat, double* eigvals,
		  rokko::parameters const& params, timer& timer) {
  timer.start(timer_id::diagonalize_diagonalize);
  //std::string matrix_part = "upper";
  char uplow = 'U';
  get_key(params, "uplow", uplow);

  int dim = mat.outerSize();
  int info;
  if(mat.is_col_major())
    info = LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'N', 'U', dim, &mat(0,0), dim, eigvals);
  else
    info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'N', 'U', dim, &mat(0,0), dim, eigvals);
  timer.stop(timer_id::diagonalize_diagonalize);
  timer.start(timer_id::diagonalize_finalize);
  if (info) {
    std::cerr << "error at dsyev function. info=" << info  << std::endl;
    exit(1);
  }
  timer.stop(timer_id::diagonalize_finalize);
  return info;
}

// dsyevd eigenvalues / eigenvectors
template<typename MATRIX_MAJOR>
int diagonalize_d(localized_matrix<double, MATRIX_MAJOR>& mat, double* eigvals,
		  localized_matrix<double, MATRIX_MAJOR>& eigvecs,
		  rokko::parameters const& params, timer& timer) {
  timer.start(timer_id::diagonalize_diagonalize);
  //std::string matrix_part = "upper";
  char uplow = 'U';
  get_key(params, "uplow", uplow);

  int dim = mat.outerSize();
  int info;
  if(mat.is_col_major())
    info = LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'U', dim, &mat(0,0), dim, eigvals);
  else
    info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', dim, &mat(0,0), dim, eigvals);
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


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_D_HPP
