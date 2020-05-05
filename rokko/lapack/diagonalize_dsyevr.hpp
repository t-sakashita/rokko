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

#ifndef ROKKO_LAPACK_DIAGONALIZE_DSYEVR_HPP
#define ROKKO_LAPACK_DIAGONALIZE_DSYEVR_HPP

#include <rokko/parameters.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>
#include <rokko/lapack/syevr.hpp>

namespace rokko {
namespace lapack {

// dsyevr only eigenvalues
template<int MATRIX_MAJOR>
parameters diagonalize_dsyevr(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, double* eigvals,
			      parameters const& params) {
  parameters params_out;
  double abstol = params.defined("abstol") ? params.get<double>("abstol") : 0.;
  params_out.set("abstol", abstol);

  lapack_int il, iu;
  double vl, vu;
  const char range = get_eigenvalues_range(params, vl, vu, il, iu);
  const char uplow = get_matrix_part(params);

  lapack_int m;  // output: found eigenvalues
  const int dim = mat.outerSize();
  std::vector<lapack_int> isuppz(2*dim+1);

  int info = syevr(range, uplow, mat,
                   vl, vu, il, iu, abstol,
                   m, eigvals, isuppz);

  if (info) {
    std::cerr << "error at dsyevr function. info=" << info << std::endl;
  }
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("isuppz", isuppz);

  if (params.get_bool("verbose")) {
    print_verbose("dsyevr", 'N', range, uplow, vl, vu, il, iu, params_out);
  }

  return params_out;
}


// dsyevr eigenvalues / eigenvectors
template<int MATRIX_MAJOR>
parameters diagonalize_dsyevr(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, double* eigvals,
			      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			      parameters const& params) {
  rokko::parameters params_out;

  double abstol = params.defined("abstol") ? params.get<double>("abstol") : 0.;
  params_out.set("abstol", abstol);

  lapack_int il = 0, iu = 0;
  double vl = 0, vu = 0;
  const char range = get_eigenvalues_range(params, vl, vu, il, iu);
  const char uplow = get_matrix_part(params);

  lapack_int m;  // output: found eigenvalues
  const int dim = mat.outerSize();
  std::vector<lapack_int> isuppz(2*dim+1);

  int info = syevr(range, uplow, mat,
                   vl, vu, il, iu, abstol,
                   m, eigvals, eigvecs, isuppz);

  if (info) {
    std::cerr << "error at dsyevr function. info=" << info << std::endl;
  }
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("isuppz", isuppz);
  
  if (params.get_bool("verbose")) {
    print_verbose("dsyevr", 'V', range, uplow, vl, vu, il, iu, params_out);
  }
  return params_out;
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_DSYEVR_HPP
