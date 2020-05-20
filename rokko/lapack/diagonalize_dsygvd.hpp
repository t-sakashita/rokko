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
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>

namespace rokko {
namespace lapack {

// dsygvd only eigenvalues
template<int MATRIX_MAJOR>
parameters diagonalize_dsygvd(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mata, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& matb,
			      double *const eigvals,
			      parameters const& params) {
  parameters params_out;
  const char jobz = 'N';  // only eigenvalues
  const char uplow = lapack::get_matrix_part(params);
  const int dim = mata.innerSize();
  const int lda = mata.outerSize();
  const int ldb = matb.outerSize();
  int info;

  if(mata.is_col_major())
    info = LAPACKE_dsygvd(LAPACK_COL_MAJOR, 1, jobz, uplow, dim, mata.data(), lda, matb.data(), ldb, eigvals);
  else
    info = LAPACKE_dsygvd(LAPACK_ROW_MAJOR, 1, jobz, uplow, dim, mata.data(), lda, matb.data(), ldb, eigvals);

  if (info) {
    std::cerr << "error at dsygvd function. info=" << info  << std::endl;
  }
  if (params.get_bool("verbose")) {
    print_verbose("sygvd", jobz, uplow);
  }

  return params_out;
}

// dsygvd eigenvalues / eigenvectors
template<int MATRIX_MAJOR>
parameters diagonalize_dsygvd(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mata, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& matb,
			      double* eigvals, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			      parameters const& params) {
  parameters params_out;
  const char jobz = 'V';  // eigenvalues / eigenvectors
  const char uplow = get_matrix_part(params);

  const int dim = mata.innerSize();
  const int lda = mata.outerSize();
  const int ldb = matb.outerSize();
  int info;

  if(mata.is_col_major())
    info = LAPACKE_dsygvd(LAPACK_COL_MAJOR, 1, jobz, uplow, dim, mata.data(), lda, matb.data(), ldb, eigvals);
  else
    info = LAPACKE_dsygvd(LAPACK_ROW_MAJOR, 1, jobz, uplow, dim, mata.data(), lda, matb.data(), ldb, eigvals);

  eigvecs = mata;
  if (info) {
    std::cerr << "error at dsygvd function. info=" << info  << std::endl;
  }
  if (params.get_bool("verbose")) {
    print_verbose("sygvd", jobz, uplow);
  }

  return params_out;
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_DSYGVD_HPP
