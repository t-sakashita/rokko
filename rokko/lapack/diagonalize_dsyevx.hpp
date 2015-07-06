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

#ifndef ROKKO_LAPACK_DIAGONALIZE_DSYEVX_HPP
#define ROKKO_LAPACK_DIAGONALIZE_DSYEVX_HPP

#include <rokko/parameters.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.h>
#include <rokko/lapack/diagonalize_get_parameters.hpp>

namespace rokko {
namespace lapack {

// dsyevx only eigenvalues
template<typename MATRIX_MAJOR>
int diagonalize_dsyevx(localized_matrix<double, MATRIX_MAJOR>& mat, double* eigvals,
		       parameters const& params, timer& timer) {
  char jobz = 'N';  // only eigenvalues
  parameters params_out;
  int dim = mat.outerSize();
  int ldim_mat = mat.innerSize();
  lapack_int m;  // output: found eigenvalues
  double abstol = 0.;  // defalut value = 0
  get_key(params, "abstol", abstol);
  params_out.set("abstol", abstol);

  lapack_int il, iu;
  double vl, vu;
  char range = get_eigenvalues_range(params, vl, vu, il, iu);
  char uplow = get_matrix_part(params);

  std::vector<lapack_int> ifail(dim);
  timer.start(timer_id::diagonalize_diagonalize);
  int info;
  if(mat.is_col_major())
    info = LAPACKE_dsyevx(LAPACK_COL_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, NULL, ldim_mat, &ifail[0]);
  else
    info = LAPACKE_dsyevx(LAPACK_ROW_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, NULL, ldim_mat, &ifail[0]);
  timer.stop(timer_id::diagonalize_diagonalize);
  timer.start(timer_id::diagonalize_finalize);
  if (info) {
    std::cerr << "error at dsyevx function. info=" << info << std::endl;
    if (info < 0) {
      std::cerr << "This means that ";
      std::cerr << "the " << abs(info) << "-th argument had an illegal value." << std::endl;
    }
    exit(1);
  }
  params_out.set("m", m);
  params_out.set("ifail", ifail);
  
  if (params.get_bool("verbose")) {
    print_verbose("dsyevx", jobz, range, uplow, vl, vu, il, iu, params_out);
  }
  timer.stop(timer_id::diagonalize_finalize);
  return info;
}


// dsyevx eigenvalues / eigenvectors
template<typename MATRIX_MAJOR>
int diagonalize_dsyevx(localized_matrix<double, MATRIX_MAJOR>& mat, double* eigvals,
		       localized_matrix<double, MATRIX_MAJOR>& eigvecs,
		       parameters const& params, timer& timer) {
  char jobz = 'V';  // eigenvalues / eigenvectors

  rokko::parameters params_out;
  int dim = mat.outerSize();
  int ldim_mat = mat.innerSize();
  int ldim_eigvec = eigvecs.innerSize();
  std::vector<lapack_int> ifail(dim);

  lapack_int m;  // output: found eigenvalues
  double abstol = 0.;  // defalut value = 0
  get_key(params, "abstol", abstol);
  params_out.set("abstol", abstol);

  lapack_int il, iu;
  double vl, vu;
  char range = get_eigenvalues_range(params, vl, vu, il, iu);
  char uplow = get_matrix_part(params);

  timer.start(timer_id::diagonalize_diagonalize);
  int info;
  if(mat.is_col_major())
    info = LAPACKE_dsyevx(LAPACK_COL_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, &eigvecs(0,0), ldim_eigvec, &ifail[0]);
  else
    info = LAPACKE_dsyevx(LAPACK_ROW_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, &eigvecs(0,0), ldim_eigvec, &ifail[0]);
  timer.stop(timer_id::diagonalize_diagonalize);
  timer.start(timer_id::diagonalize_finalize);
  if (info) {
    std::cerr << "Error at dsyevx function. info=" << info << std::endl;
    if (params.get_bool("verbose")) {
      std::cerr << "This means that ";
      if (info < 0) {
	std::cerr << "the " << abs(info) << "-th argument had an illegal value." << std::endl;
      } else {
	std::cerr << "This means that "	<< info << " eigenvectors failed to converge." << std::endl;
	std::cerr << "The indices of the eigenvectors that failed to converge:" << std::endl;
	for (int i=0; i<ifail.size(); ++i) {
	  if (ifail[i] == 0) break;
	  std::cerr << ifail[i] << " ";
	}
	std::cerr << std::endl;
      }
    }
    exit(1);
  }
  params_out.set("m", m);
  params_out.set("ifail", ifail);
  
  if (params.get_bool("verbose")) {
    print_verbose("dsyevx", jobz, range, uplow, vl, vu, il, iu, params_out);
  }
  timer.stop(timer_id::diagonalize_finalize);
  return info;
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_DSYEVX_HPP
