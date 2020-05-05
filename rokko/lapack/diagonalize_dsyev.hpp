/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_LAPACK_DIAGONALIZE_DSYEV_HPP
#define ROKKO_LAPACK_DIAGONALIZE_DSYEV_HPP

#include <rokko/parameters.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>
#include <rokko/lapack/syev.hpp>

namespace rokko {
namespace lapack {

// dsyev only eigenvalues
template<int MATRIX_MAJOR>
parameters diagonalize_dsyev(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, double* eigvals,
			     parameters const& params) {
  rokko::parameters params_out;
  const char jobz = 'N';  // only eigenvalues
  char uplow = get_matrix_part(params);

  int info = syev(jobz, uplow, mat, eigvals);

  params_out.set("info", info);
  if (info) {
    std::cerr << "error at dsyev function. info=" << info  << std::endl;
  }
  if (params.get_bool("verbose")) {
    print_verbose("dsyev", jobz, uplow);
  }

  return params_out;
}

// dsyev eigenvalues / eigenvectors
template<int MATRIX_MAJOR>
parameters diagonalize_dsyev(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, double* eigvals,
			     Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			     parameters const& params) {
  rokko::parameters params_out;
  const char jobz = 'V';  // eigenvalues / eigenvectors
  char uplow = get_matrix_part(params);

  int info = syev(jobz, uplow, mat, eigvals);
  eigvecs = mat;

  params_out.set("info", info);
  if (info) {
    std::cerr << "error at dsyev function. info=" << info  << std::endl;
  }
  if (params.get_bool("verbose")) {
    print_verbose("dsyev", jobz, uplow);
  }

  return params_out;
}

} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_DSYEV_HPP
