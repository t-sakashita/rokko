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

#ifndef ROKKO_LAPACK_DIAGONALIZE_DSYGVD_HPP
#define ROKKO_LAPACK_DIAGONALIZE_DSYGVD_HPP

#include <rokko/parameters.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>

namespace rokko {
namespace lapack {

// dsygvd only eigenvalues
template<typename MATRIX_MAJOR>
parameters diagonalize_dsygvd(localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
			      double* eigvals,
			      parameters const& params) {
  parameters params_out;
  char jobz = 'N';  // only eigenvalues
  char uplow = lapack::get_matrix_part(params);
  int dim = mata.innerSize();
  int lda = mata.outerSize();
  int ldb = matb.outerSize();
  int info;

  if(mata.is_col_major())
    info = LAPACKE_dsygvd(LAPACK_COL_MAJOR, 1, jobz, uplow, dim, &mata(0,0), lda, &matb(0,0), ldb, eigvals);
  else
    info = LAPACKE_dsygvd(LAPACK_ROW_MAJOR, 1, jobz, uplow, dim, &mata(0,0), lda, &matb(0,0), ldb, eigvals);

  if (info) {
    std::cerr << "error at dsygvd function. info=" << info  << std::endl;
    exit(1);
  }
  if (params.get_bool("verbose")) {
    print_verbose("dsygvd", jobz, uplow);
  }

  return params_out;
}

// dsygvd eigenvalues / eigenvectors
template<typename MATRIX_MAJOR>
parameters diagonalize_dsygvd(localized_matrix<double, MATRIX_MAJOR>& mata, localized_matrix<double, MATRIX_MAJOR>& matb,
			      double* eigvals, localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			      parameters const& params) {
  parameters params_out;
  char jobz = 'V';  // eigenvalues / eigenvectors
  char uplow = get_matrix_part(params);

  int dim = mata.innerSize();
  int lda = mata.outerSize();
  int ldb = matb.outerSize();
  int info;

  if(mata.is_col_major())
    info = LAPACKE_dsygvd(LAPACK_COL_MAJOR, 1, jobz, uplow, dim, &mata(0,0), lda, &matb(0,0), ldb, eigvals);
  else
    info = LAPACKE_dsygvd(LAPACK_ROW_MAJOR, 1, jobz, uplow, dim, &mata(0,0), lda, &matb(0,0), ldb, eigvals);

  eigvecs = mata;
  if (info) {
    std::cerr << "error at dsygvd function. info=" << info  << std::endl;
    exit(1);
  }
  if (params.get_bool("verbose")) {
    print_verbose("dsygvd", jobz, uplow);
  }

  return params_out;
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_DSYGVD_HPP
