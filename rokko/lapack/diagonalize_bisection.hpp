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
parameters diagonalize_bisection(localized_matrix<double, MATRIX_MAJOR>& mat, double* eigvals,
				 rokko::parameters const& params) {
  rokko::parameters params_out;
  char jobz = 'N';  // only eigenvalues
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
  params_out.set("abstol", abstol);
  char uplow = get_matrix_part(params);

  lapack_int il, iu;
  double vl, vu;
  char range = get_eigenvalues_range(params, vl, vu, il, iu);

  std::vector<lapack_int> ifail(dim);
  int info;
  if(mat.is_col_major())
    info = LAPACKE_dsyevx(LAPACK_COL_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, NULL, ldim_mat, &ifail[0]);
  else
    info = LAPACKE_dsyevx(LAPACK_ROW_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, NULL, ldim_mat, &ifail[0]);
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
    print_verbose("dsyevx (bisection)", jobz, range, uplow, vl, vu, il, iu, params_out);
  }

  return params_out;
}


// dsyevx eigenvalues / eigenvectors
template<typename MATRIX_MAJOR>
parameters diagonalize_bisection(localized_matrix<double, MATRIX_MAJOR>& mat, double* eigvals,
				 localized_matrix<double, MATRIX_MAJOR>& eigvecs,
				 parameters const& params) {
  rokko::parameters params_out;
  char jobz = 'V';  // eigenvalues / eigenvectors
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
  params_out.set("abstol", abstol);
  char uplow = get_matrix_part(params);

  lapack_int il, iu;
  double vl, vu;
  char range = get_eigenvalues_range(params, vl, vu, il, iu);

  int info;
  if(mat.is_col_major())
    info = LAPACKE_dsyevx(LAPACK_COL_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, &eigvecs(0,0), ldim_eigvec, &ifail[0]);
  else
    info = LAPACKE_dsyevx(LAPACK_ROW_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, &eigvecs(0,0), ldim_eigvec, &ifail[0]);

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
    print_verbose("dsyevx (bisecition)", jobz, range, uplow, vl, vu, il, iu, params_out);
  }

  return params_out;
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_BISECTION_HPP
