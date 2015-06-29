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

#ifndef ROKKO_LAPACK_DIAGONALIZE_HPP
#define ROKKO_LAPACK_DIAGONALIZE_HPP

#include <rokko/parameters.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.h>
#include <rokko/lapack/diagonalize_get_parameters.hpp>

namespace rokko {
namespace lapack {

// dsyev only eigenvalues
template<typename MATRIX_MAJOR>
int diagonalize(localized_matrix<double, MATRIX_MAJOR>& mat, double* eigvals,
		rokko::parameters const& params, timer& timer) {
  char jobz = 'N';  // only eigenvalues

  std::string matrix_part = "upper"; // default is "upper"
  char uplow = 'U';
  get_matrix_part(params, matrix_part, uplow);

  timer.start(timer_id::diagonalize_diagonalize);
  int dim = mat.outerSize();
  int info;
  if(mat.is_col_major())
    info = LAPACKE_dsyev(LAPACK_COL_MAJOR, jobz, uplow, dim, &mat(0,0), dim, eigvals);
  else
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, uplow, dim, &mat(0,0), dim, eigvals);
  timer.stop(timer_id::diagonalize_diagonalize);
  timer.start(timer_id::diagonalize_finalize);
  if (info) {
    std::cerr << "error at dsyev function. info=" << info  << std::endl;
    exit(1);
  }
  if (params.get_bool("verbose")) {
    std::cout << "All eigenvalues/eigenvectors were requested by dsyev" << std::endl;
  }
  timer.stop(timer_id::diagonalize_finalize);
  return info;
}

// dsyev eigenvalues / eigenvectors
template<typename MATRIX_MAJOR>
int diagonalize(localized_matrix<double, MATRIX_MAJOR>& mat, double* eigvals,
                localized_matrix<double, MATRIX_MAJOR>& eigvecs,
		rokko::parameters const& params, timer& timer) {
  char jobz = 'V';  // eigenvalues / eigenvectors

  std::string matrix_part = "upper"; // default is "upper"
  char uplow = 'U';
  get_matrix_part(params, matrix_part, uplow);

  timer.start(timer_id::diagonalize_diagonalize);
  int dim = mat.outerSize();
  int info;
  if(mat.is_col_major())
    info = LAPACKE_dsyev(LAPACK_COL_MAJOR, jobz, uplow, dim, &mat(0,0), dim, eigvals);
  else
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, uplow, dim, &mat(0,0), dim, eigvals);
  timer.stop(timer_id::diagonalize_diagonalize);
  timer.start(timer_id::diagonalize_finalize);
  eigvecs = mat;
  if (info) {
    std::cerr << "error at dsyev function. info=" << info  << std::endl;
    exit(1);
  }
  if (params.get_bool("verbose")) {
    std::cout << "All eigenvalues/eigenvectors were requested by dsyev" << std::endl;
  }
  timer.stop(timer_id::diagonalize_finalize);
  return info;
}

} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_HPP
