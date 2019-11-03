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

#ifndef ROKKO_LAPACK_DIAGONALIZE_DSYEVR_HPP
#define ROKKO_LAPACK_DIAGONALIZE_DSYEVR_HPP

#include <rokko/parameters.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>

namespace rokko {
namespace lapack {

// dsyevr only eigenvalues
template<int MATRIX_MAJOR>
parameters diagonalize_dsyevr(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, double* eigvals,
			      parameters const& params) {
  parameters params_out;
  const char jobz = 'N';  // only eigenvalues
  const int dim = mat.outerSize();
  const int ld_mat = mat.innerSize();
  double abstol = 0.;  // defalut value = 0
  get_key(params, "abstol", abstol);
  params_out.set("abstol", abstol);

  lapack_int il, iu;
  double vl, vu;
  const char range = get_eigenvalues_range(params, vl, vu, il, iu);
  const char uplow = get_matrix_part(params);

  lapack_int m;  // output: found eigenvalues
  std::vector<lapack_int> isuppz(2*dim+1);
  int info;
  if(MATRIX_MAJOR == Eigen::ColMajor)
    info = LAPACKE_dsyevr(LAPACK_COL_MAJOR, jobz, range, uplow, dim,
			  &mat(0,0), ld_mat, vl, vu, il, iu,
			  abstol, &m, eigvals, NULL, ld_mat, &isuppz[0]);
  else
    info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, jobz, range, uplow, dim,
			  &mat(0,0), ld_mat, vl, vu, il, iu,
			  abstol, &m, eigvals, NULL, ld_mat, &isuppz[0]);

  if (info) {
    std::cerr << "error at dsyevr function. info=" << info << std::endl;
  }
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("isuppz", isuppz);

  if (params.get_bool("verbose")) {
    print_verbose("dsyevr", jobz, range, uplow, vl, vu, il, iu, params_out);
  }

  return params_out;
}


// dsyevr eigenvalues / eigenvectors
template<int MATRIX_MAJOR>
parameters diagonalize_dsyevr(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, double* eigvals,
			      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			      parameters const& params) {
  rokko::parameters params_out;
  const char jobz = 'V';  // eigenvalues / eigenvectors

  const int dim = mat.outerSize();
  const int ld_mat = mat.innerSize();
  const int ld_eigvecs = eigvecs.innerSize();

  double abstol = 0.;  // defalut value = 0
  get_key(params, "abstol", abstol);
  params_out.set("abstol", abstol);

  lapack_int il = 0, iu = 0;
  double vl = 0, vu = 0;
  const char range = get_eigenvalues_range(params, vl, vu, il, iu);
  const char uplow = get_matrix_part(params);

  lapack_int m;  // output: found eigenvalues
  std::vector<lapack_int> isuppz(2*dim+1);

  int info;
  if(MATRIX_MAJOR == Eigen::ColMajor)
    info = LAPACKE_dsyevr(LAPACK_COL_MAJOR, jobz, range, uplow, dim,
			  &mat(0,0), ld_mat, vl, vu, il, iu,
			  abstol, &m, eigvals, &eigvecs(0,0), ld_eigvecs, &isuppz[0]);
  else
    info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, jobz, range, uplow, dim,
			  &mat(0,0), ld_mat, vl, vu, il, iu,
			  abstol, &m, eigvals, &eigvecs(0,0), ld_eigvecs, &isuppz[0]);

  if (info) {
    std::cerr << "error at dsyevr function. info=" << info << std::endl;
  }
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("isuppz", isuppz);
  
  if (params.get_bool("verbose")) {
    print_verbose("dsyevr", jobz, range, uplow, vl, vu, il, iu, params_out);
  }
  return params_out;
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_DSYEVR_HPP
