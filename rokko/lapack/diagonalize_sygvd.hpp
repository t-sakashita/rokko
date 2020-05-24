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

#ifndef ROKKO_LAPACK_DIAGONALIZE_SYGVD_HPP
#define ROKKO_LAPACK_DIAGONALIZE_SYGVD_HPP

#include <rokko/parameters.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>

namespace rokko {
namespace lapack {

// only eigenvalues
template<typename T, int MATRIX_MAJOR, typename VEC>
parameters diagonalize_sygvd(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mata, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& matb,
			      VEC& eigvals,
			      parameters const& params) {
  parameters params_out;
  const char jobz = 'N';  // only eigenvalues
  const char uplow = lapack::get_matrix_part(params);

  constexpr int itype = 1;
  int info = sygvd(itype, jobz, uplow, mata, matb, eigvals);

  if (info) {
    std::cerr << "error at sygvd function. info=" << info << std::endl;
  }
  if (params.get_bool("verbose")) {
    print_verbose("sygvd", jobz, uplow);
  }

  return params_out;
}

// eigenvalues / eigenvectors
template<typename T, int MATRIX_MAJOR, typename VEC>
parameters diagonalize_sygvd(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mata, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& matb,
			      VEC& eigvals, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			      parameters const& params) {
  parameters params_out;
  const char jobz = 'V';  // eigenvalues / eigenvectors
  const char uplow = get_matrix_part(params);

  constexpr int itype = 1;
  int info = sygvd(itype, jobz, uplow, mata, matb, eigvals);
  eigvecs = mata;
  if (info) {
    std::cerr << "error at sygvd function. info=" << info << std::endl;
  }
  if (params.get_bool("verbose")) {
    print_verbose("sygvd", jobz, uplow);
  }

  return params_out;
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_SYGVD_HPP
