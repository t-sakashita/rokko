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

#ifndef ROKKO_LAPACK_DIAGONALIZE_DSYEVD_HPP
#define ROKKO_LAPACK_DIAGONALIZE_DSYEVD_HPP

#include <rokko/parameters.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>

namespace rokko {
namespace lapack {

// dsyevd only eigenvalues
template<int MATRIX_MAJOR>
parameters diagonalize_dsyevd(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, double* eigvals,
			      parameters const& params) {
  rokko::parameters params_out;
  const char jobz = 'N';  // only eigenvalues
  const char uplow = lapack::get_matrix_part(params);

  const int dim = mat.outerSize();
  const int ldim = mat.innerSize();
  int info;

  if(MATRIX_MAJOR == Eigen::ColMajor)
    info = LAPACKE_dsyevd(LAPACK_COL_MAJOR, jobz, uplow, dim, &mat(0,0), ldim, eigvals);
  else
    info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, jobz, uplow, dim, &mat(0,0), ldim, eigvals);

  params_out.set("info", info);
  if (info) {
    std::cerr << "error at dsyevd function. info=" << info  << std::endl;
    //exit(1);
  }
  if (params.get_bool("verbose")) {
    print_verbose("dsyevd", jobz, uplow);
  }

  return params_out;
}

// dsyevd eigenvalues / eigenvectors
template<int MATRIX_MAJOR>
parameters diagonalize_dsyevd(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, double* eigvals,
			      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			      parameters const& params) {
  rokko::parameters params_out;
  const char jobz = 'V';  // eigenvalues / eigenvectors
  const char uplow = get_matrix_part(params);
  const int dim = mat.outerSize();
  const int ldim = mat.innerSize();
  int info;

  if(MATRIX_MAJOR == Eigen::ColMajor)
    info = LAPACKE_dsyevd(LAPACK_COL_MAJOR, jobz, uplow, dim, &mat(0,0), ldim, eigvals);
  else
    info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, jobz, uplow, dim, &mat(0,0), ldim, eigvals);

  eigvecs = mat;
  params_out.set("info", info);
  if (info) {
    std::cerr << "error at dsyevd function. info=" << info  << std::endl;
    exit(1);
  }
  if (params.get_bool("verbose")) {
    print_verbose("dsyevd", jobz, uplow);
  }

  return params_out;
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_DSYEVD_HPP
