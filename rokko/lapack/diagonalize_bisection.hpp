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

#ifndef ROKKO_LAPACK_DIAGONALIZE_BISECTION_HPP
#define ROKKO_LAPACK_DIAGONALIZE_BISECTION_HPP

#include <rokko/parameters.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.h>

namespace rokko {
namespace lapack {

// dsyevx only eigenvalues
template<typename MATRIX_MAJOR>
int diagonalize_bisection(localized_matrix<double, MATRIX_MAJOR>& mat, double* eigvals,
			  rokko::parameters const& params, timer& timer) {
  char jobz = 'N';  // only eigenvalues
  rokko::parameters params_out;
  int dim = mat.outerSize();
  int ldim_mat = mat.innerSize();
  lapack_int m;  // output: found eigenvalues
  double abstol;
  get_key(params, "abstol", abstol);
  if (abstol < 0) {
    std::cerr << "Error in diagonalize_bisection" << std::endl
	      << "abstol is negative value, which means QR method." << std::endl
	      << "To use dsyevx as bisection solver, set abstol a positive value" << std::endl;
    throw;
  }
  if (!params.defined("abstol")) {  // default: optimal value for bisection method
    abstol = 2 * LAPACKE_dlamch('S');
  }

  std::string matrix_part = "upper"; // default is "upper"
  char uplow = 'U';
  get_matrix_part(params, matrix_part, uplow);

  char range = 'A';  // default is 'A'
  lapack_int il = 0, iu = 0;
  double vl = 0, vu = 0;
  bool is_upper_value, is_upper_index, is_lower_value, is_lower_index;
  get_eigenvalues_range(params, matrix_part, range, vu, vl, iu, il, is_upper_value, is_lower_value, is_upper_index, is_lower_index);

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
    if (range == 'A')
      std::cout << "All eigenvalues were requested by dsyevx (bisection)" << std::endl;
    else if (is_upper_value && is_lower_value)	    
      std::cout << "Eigenvalues/eigenvectors contained in the interval [" << vl << ", " << vu << "]" << " were requested" << std::endl;
    else if (is_upper_index && is_lower_index)
      std::cout << "Eigenvalues/eigenvectors from " << il << "th" << " to " << iu << "th" << " were requested" << std::endl;
    std::cout << "The number of found eigenvalues are " << m << std::endl;
    std::cout << "The " << matrix_part << " part of the matrix was used" << std::endl;
    if (!params.defined("abstol")) {
      std::cout << "abstol was not specified, so used optimal value for bisection method: 2 * LAPACKE_dlamch('S')" << std::endl;
    }
    std::cout << "abstol=" << abstol << std::endl;
  }
  timer.stop(timer_id::diagonalize_finalize);
  return info;
}


// dsyevx eigenvalues / eigenvectors
template<typename MATRIX_MAJOR>
int diagonalize_bisection(localized_matrix<double, MATRIX_MAJOR>& mat, double* eigvals,
			  localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			  rokko::parameters const& params, timer& timer) {
  char jobz = 'V';  // eigenvalues / eigenvectors

  rokko::parameters params_out;
  int dim = mat.outerSize();
  int ldim_mat = mat.innerSize();
  int ldim_eigvec = eigvecs.innerSize();
  std::vector<lapack_int> ifail(dim);

  lapack_int m;  // output: found eigenvalues
  double abstol;
  get_key(params, "abstol", abstol);
  if (abstol < 0) {
    std::cerr << "Error in diagonalize_bisection" << std::endl
	      << "abstol is negative value, which means QR method." << std::endl
	      << "To use dsyevx as bisection solver, set abstol a positive value" << std::endl;
    throw;
  }
  if (!params.defined("abstol")) {  // default: optimal value for bisection method
    abstol = 2 * LAPACKE_dlamch('S');
  }

  std::string matrix_part = "upper"; // default is "upper"
  char uplow = 'U';
  get_matrix_part(params, matrix_part, uplow);

  char range = 'A';  // default is 'A'
  lapack_int il = 0, iu = 0;
  double vl = 0, vu = 0;
  bool is_upper_value, is_upper_index, is_lower_value, is_lower_index;
  get_eigenvalues_range(params, matrix_part, range, vu, vl, iu, il, is_upper_value, is_lower_value, is_upper_index, is_lower_index);

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
    if (range == 'A')
      std::cout << "All eigenvalues/eigenvectors were requested by dsyevx (bisection)" << std::endl;
    else if (is_upper_value && is_lower_value)
      std::cout << "Eigenvalues/eigenvectors contained in the interval [" << vl << ", " << vu << "]" << " were requested" << std::endl;
    else if (is_upper_index && is_lower_index)
      std::cout << "Eigenvalues/eigenvectors from " << il << "th" << " to " << iu << "th" << " were requested" << std::endl;
    std::cout << "The number of found eigenvalues are " << m << std::endl;
    std::cout << "The " << matrix_part << " part of the matrix was used" << std::endl;
    if (!params.defined("abstol")) {
      std::cout << "abstol was not specified, so used optimal value for bisection method: 2 * LAPACKE_dlamch('S')" << std::endl;
    }
    std::cout << "abstol=" << abstol << std::endl;
  }
  timer.stop(timer_id::diagonalize_finalize);
  return info;
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_BISECTION_HPP
