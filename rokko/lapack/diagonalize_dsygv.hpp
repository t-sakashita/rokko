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

#ifndef ROKKO_LAPACK_DIAGONALIZE_DSYGV_HPP
#define ROKKO_LAPACK_DIAGONALIZE_DSYGV_HPP

#include <rokko/parameters.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>

namespace rokko {
namespace lapack {

// dsygv only eigenvalues
template<int MATRIX_MAJOR>
parameters diagonalize_dsygv(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mata, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& matb,
			     double* eigvals,
			     parameters const& params) {
  parameters params_out;
  char jobz = 'N';  // only eigenvalues
  char uplow = get_matrix_part(params);

  int dim = mata.innerSize();
  int lda = mata.outerSize();
  int ldb = matb.outerSize();
  int info;

  if(mata.is_col_major())
    info = LAPACKE_dsygv(LAPACK_COL_MAJOR, 1, jobz, uplow, dim, &mata(0,0), lda, &matb(0,0), ldb, eigvals);
  else
    info = LAPACKE_dsygv(LAPACK_ROW_MAJOR, 1, jobz, uplow, dim, &mata(0,0), lda, &matb(0,0), ldb, eigvals);

  if (info) {
    std::cerr << "error at dsygv function. info=" << info  << std::endl;
    exit(1);
  }
  if (params.get_bool("verbose")) {
    print_verbose("dsygv", jobz, uplow);
  }

  return params_out;
}

// dsygv eigenvalues / eigenvectors
template<int MATRIX_MAJOR>
parameters diagonalize_dsygv(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mata, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& matb,
			     double* eigvals,
			     Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			     parameters const& params) {
  parameters params_out;
  char jobz = 'V';  // eigenvalues / eigenvectors
  char uplow = get_matrix_part(params);

  int dim = mata.innerSize();
  int lda = mata.outerSize();
  int ldb = matb.outerSize();
  int info;

  if(mata.is_col_major())
    info = LAPACKE_dsygv(LAPACK_COL_MAJOR, 1, jobz, uplow, dim, &mata(0,0), lda, &matb(0,0), ldb, eigvals);
  else
    info = LAPACKE_dsygv(LAPACK_ROW_MAJOR, 1, jobz, uplow, dim, &mata(0,0), lda, &matb(0,0), ldb, eigvals);

  eigvecs = mata;
  if (info) {
    std::cerr << "error at dsygv function. info=" << info  << std::endl;
    exit(1);
  }
  if (params.get_bool("verbose")) {
    print_verbose("dsygv", jobz, uplow);
  }

  return params_out;
}

} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_DSYGV_HPP
