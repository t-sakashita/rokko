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

#pragma once

#include <rokko/parameters.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>
#include <rokko/lapack/syev.hpp>

namespace rokko {
namespace lapack {

// only eigenvalues
template<typename T, int MATRIX_MAJOR, typename VEC>
parameters diagonalize_syev(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, VEC& eigvals,
			     parameters const& params) {
  rokko::parameters params_out;
  const char jobz = 'N';  // only eigenvalues
  const char uplow = get_matrix_part(params);

  const auto info = syev(jobz, uplow, mat, eigvals);

  params_out.set("info", info);
  if (info) {
    std::cerr << "error at syev function. info=" << info  << std::endl;
  }
  if (params.get_bool("verbose")) {
    print_verbose("syev", jobz, uplow);
  }

  return params_out;
}

// eigenvalues / eigenvectors
template<typename T, int MATRIX_MAJOR, typename VEC>
parameters diagonalize_syev(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, VEC& eigvals,
			     Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			     parameters const& params) {
  rokko::parameters params_out;
  const char jobz = 'V';  // eigenvalues / eigenvectors
  const char uplow = get_matrix_part(params);

  const auto info = syev(jobz, uplow, mat, eigvals);
  eigvecs = mat;

  params_out.set("info", info);
  if (info) {
    std::cerr << "error at syev function. info=" << info  << std::endl;
  }
  if (params.get_bool("verbose")) {
    print_verbose("syev", jobz, uplow);
  }

  return params_out;
}

} // namespace lapack
} // namespace rokko
