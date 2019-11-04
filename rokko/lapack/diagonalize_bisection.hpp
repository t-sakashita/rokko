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
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>

namespace rokko {
namespace lapack {

// dsyevx only eigenvalues
template<int MATRIX_MAJOR>
parameters diagonalize_bisection(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, double* eigvals,
				 rokko::parameters const& params) {
  rokko::parameters params_out;
  const char jobz = 'N';  // only eigenvalues
  const int dim = mat.outerSize();
  const int ldim_mat = mat.innerSize();
  lapack_int m;  // output: found eigenvalues
  double abstol = params.defined("abstol") ? params.get<double>("abstol") : 2*LAPACKE_dlamch('S');
  if (abstol < 0) {
    std::stringstream msg;
    msg << "lapack::diagonalize_bisection() : " << std::endl
	<< "abstol is negative value, which means QR method." << std::endl
	<< "To use dsyevx as bisection solver, set abstol a positive value" << std::endl;
    throw std::invalid_argument(msg.str());
  }
  params_out.set("abstol", abstol);
  const char uplow = get_matrix_part(params);

  lapack_int il, iu;
  double vl, vu;
  const char range = get_eigenvalues_range(params, vl, vu, il, iu);

  std::vector<lapack_int> ifail(dim);
  int info;
  if(MATRIX_MAJOR == Eigen::ColMajor)
    info = LAPACKE_dsyevx(LAPACK_COL_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, NULL, ldim_mat, &ifail[0]);
  else
    info = LAPACKE_dsyevx(LAPACK_ROW_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, NULL, ldim_mat, &ifail[0]);
  if (info) {
    std::stringstream msg;
    msg << "lapack::diagonalize_bisection() : " << std::endl
	<< "error at dsyevx function. info=" << info << std::endl;
    if (info < 0) {
      msg << "This means that "
	  << "the " << abs(info) << "-th argument had an illegal value." << std::endl;
    }
    throw std::invalid_argument(msg.str());
  }
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("ifail", ifail);
  
  if (params.get_bool("verbose")) {
    print_verbose("dsyevx (bisection)", jobz, range, uplow, vl, vu, il, iu, params_out);
  }

  return params_out;
}


// dsyevx eigenvalues / eigenvectors
template<int MATRIX_MAJOR>
parameters diagonalize_bisection(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, double* eigvals,
				 Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
				 parameters const& params) {
  rokko::parameters params_out;
  const char jobz = 'V';  // eigenvalues / eigenvectors
  const int dim = mat.outerSize();
  const int ldim_mat = mat.innerSize();
  const int ldim_eigvec = eigvecs.innerSize();
  std::vector<lapack_int> ifail(dim);

  lapack_int m;  // output: found eigenvalues
  double abstol = params.defined("abstol") ? params.get<double>("abstol") : 2*LAPACKE_dlamch('S');
  if (abstol < 0) {
    std::stringstream msg;
    msg << "lapack::diagonalize_bisection() : " << std::endl
	<< "abstol is negative value, which means QR method." << std::endl
	<< "To use dsyevx as bisection solver, set abstol a positive value" << std::endl;
    throw std::invalid_argument(msg.str());
  }
  params_out.set("abstol", abstol);
  char uplow = get_matrix_part(params);

  lapack_int il, iu;
  double vl, vu;
  const char range = get_eigenvalues_range(params, vl, vu, il, iu);

  int info;
  if(MATRIX_MAJOR == Eigen::ColMajor)
    info = LAPACKE_dsyevx(LAPACK_COL_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, &eigvecs(0,0), ldim_eigvec, &ifail[0]);
  else
    info = LAPACKE_dsyevx(LAPACK_ROW_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, &eigvecs(0,0), ldim_eigvec, &ifail[0]);

  if (info) {
    std::stringstream msg;
    msg << "lapack::diagonalize_bisection() : "
	<< "error at dsyevx function. info=" << info << std::endl;
    if (params.get_bool("verbose")) {
      msg << "This means that ";
      if (info < 0) {
	msg << "the " << abs(info) << "-th argument had an illegal value." << std::endl;
      } else {
	msg << "This means that " << info << " eigenvectors failed to converge." << std::endl;
	msg << "The indices of the eigenvectors that failed to converge:" << std::endl;
	for (std::size_t i = 0; i < ifail.size(); ++i) {
	  if (ifail[i] == 0) break;
	  msg << ifail[i] << " ";
	}
	msg << std::endl;
      }
    }
    throw std::invalid_argument(msg.str());
  }
  params_out.set("info", info);
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
