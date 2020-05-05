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

#ifndef ROKKO_LAPACK_DIAGONALIZE_SYEVD_HPP
#define ROKKO_LAPACK_DIAGONALIZE_SYEVD_HPP

#include <rokko/parameters.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>
#include <rokko/lapack/syevd.hpp>

namespace rokko {
namespace lapack {

// only eigenvalues
template<typename T, int MATRIX_MAJOR>
parameters diagonalize_syevd(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, T* eigvals,
			      parameters const& params) {
  rokko::parameters params_out;
  const char uplow = lapack::get_matrix_part(params);

  int info = syevd('N', uplow, mat, eigvals);

  params_out.set("info", info);
  if (info) {
    std::cerr << "error at syevd function. info=" << info  << std::endl;
  }
  if (params.get_bool("verbose")) {
    print_verbose("syevd", 'N', uplow);
  }

  return params_out;
}

// eigenvalues / eigenvectors
template<typename T, int MATRIX_MAJOR>
parameters diagonalize_syevd(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, T* eigvals,
			      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			      parameters const& params) {
  rokko::parameters params_out;
  const char uplow = get_matrix_part(params);

  int info = syevd('V', uplow, mat, eigvals);
  eigvecs = mat;

  params_out.set("info", info);
  if (info) {
    std::cerr << "error at syevd function. info=" << info  << std::endl;
  }
  if (params.get_bool("verbose")) {
    print_verbose("syevd", 'V', uplow);
  }

  return params_out;
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_SYEVD_HPP
