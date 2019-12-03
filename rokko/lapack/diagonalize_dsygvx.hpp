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

#ifndef ROKKO_LAPACK_DIAGONALIZE_DSYGVX_HPP
#define ROKKO_LAPACK_DIAGONALIZE_DSYGVX_HPP

#include <rokko/parameters.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>

namespace rokko {
namespace lapack {

// dsygvx only eigenvalues
template<int MATRIX_MAJOR>
parameters diagonalize_dsygvx(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mata, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& matb,
			      double* eigvals,
			      parameters const& params) {
  parameters params_out;
  const char jobz = 'N';  // only eigenvalues
  const int dim = mata.innerSize();
  const int lda = mata.outerSize();
  const int ldb = matb.outerSize();
  lapack_int m;  // output: found eigenvalues
  double abstol = params.defined("abstol") ? params.get<double>("abstol") : 0.;
  params_out.set("abstol", abstol);

  lapack_int il, iu;
  double vl, vu;
  const char range = get_eigenvalues_range(params, vl, vu, il, iu);
  const char uplow = get_matrix_part(params);

  std::vector<lapack_int> ifail(dim);
  int info;
  if(mata.is_col_major())
    info = LAPACKE_dsygvx(LAPACK_COL_MAJOR, 1, jobz, range, uplow, dim,
			  mata.data(), lda, matb.data(), ldb, vl, vu, il, iu,
			  abstol, &m, eigvals, NULL, lda, ifail.data());
  else
    info = LAPACKE_dsygvx(LAPACK_ROW_MAJOR, 1, jobz, range, uplow, dim,
			  mata.data(), lda, matb.data(), ldb, vl, vu, il, iu,
			  abstol, &m, eigvals, NULL, lda, ifail.data());

  if (info) {
    std::cerr << "error at dsygvx function. info=" << info << std::endl;
    if (info < 0) {
      std::cerr << "This means that ";
      std::cerr << "the " << abs(info) << "-th argument had an illegal value." << std::endl;
    }
  }
  params_out.set("m", m);
  params_out.set("ifail", ifail);
  
  if (params.get_bool("verbose")) {
    print_verbose("dsygvx", jobz, range, uplow, vl, vu, il, iu, params_out);
  }

  return params_out;
}


// dsygvx eigenvalues / eigenvectors
template<int MATRIX_MAJOR>
parameters diagonalize_dsygvx(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mata, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& matb,
			      double* eigvals,
			      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			      parameters const& params) {
  rokko::parameters params_out;
  const char jobz = 'V';  // eigenvalues / eigenvectors
  const int dim = mata.innerSize();
  const int lda = mata.outerSize();
  const int ldb = matb.outerSize();
  const int ldim_eigvec = eigvecs.innerSize();
  std::vector<lapack_int> ifail(dim);

  lapack_int m;  // output: found eigenvalues
  double abstol = params.defined("abstol") ? params.get<double>("abstol") : 0.;
  params_out.set("abstol", abstol);

  lapack_int il, iu;
  double vl, vu;
  const char range = get_eigenvalues_range(params, vl, vu, il, iu);
  const char uplow = get_matrix_part(params);

  int info;
  if(mata.is_col_major())
    info = LAPACKE_dsygvx(LAPACK_COL_MAJOR, 1, jobz, range, uplow, dim,
			  mata.data(), lda, matb.data(), ldb, vl, vu, il, iu,
			  abstol, &m, eigvals, eigvecs.data(), ldim_eigvec, ifail.data());
  else
    info = LAPACKE_dsygvx(LAPACK_ROW_MAJOR, 1, jobz, range, uplow, dim,
			  mata.data(), lda, matb.data(), ldb, vl, vu, il, iu,
			  abstol, &m, eigvals, eigvecs.data(), ldim_eigvec, ifail.data());

  if (info) {
    std::cerr << "Error at dsygvx function. info=" << info << std::endl;
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
  }
  params_out.set("m", m);
  params_out.set("ifail", ifail);
  
  if (params.get_bool("verbose")) {
    print_verbose("dsygvx", jobz, range, uplow, vl, vu, il, iu, params_out);
  }

  return params_out;
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_DSYGVX_HPP
