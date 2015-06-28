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

#ifndef ROKKO_LAPACK_DIAGONALIZE_R_HPP
#define ROKKO_LAPACK_DIAGONALIZE_R_HPP

#include <rokko/parameters.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.h>

namespace rokko {
namespace lapack {

// dsyevr only eigenvalues
template<typename MATRIX_MAJOR>
int diagonalize_r(localized_matrix<double, MATRIX_MAJOR>& mat, double* eigvals,
		  rokko::parameters const& params, timer& timer) {
  char jobz = 'N';  // only eigenvalues
  rokko::parameters params_out;
  int dim = mat.outerSize();
  int ldim_mat = mat.innerSize();
  lapack_int m;  // output: found eigenvalues
  double abstol = 0.;  // defalut value = 0
  get_key(params, "abstol", abstol);

  std::string matrix_part = "upper"; // default is "upper"
  char uplow = 'U';
  if (params.defined("uplow"))
    matrix_part = params.get_string("uplow");
  if (params.defined("matrix_part"))
    matrix_part = params.get_string("matrix_part");
  if ((matrix_part == "upper") || (matrix_part == "U"))
    matrix_part = "upper";  uplow = 'U';
  if ((matrix_part == "lower") || (matrix_part == "L"))
    matrix_part = "lower";  uplow = 'L';

  char range = 'A';  // default is 'A'
  lapack_int il = 0, iu = 0;
  double vl = 0, vu = 0;
  bool upper_limit_double = get_key(params, "upper_limit", vu);
  bool upper_limit_int = get_key(params, "upper_limit", iu);
  bool lower_limit_double = get_key(params, "lower_limit", vl);
  bool lower_limit_int = get_key(params, "lower_limit", il);
  if (upper_limit_int && lower_limit_int)   range = 'I';
  if (upper_limit_double && lower_limit_double)   range = 'V';
  if (upper_limit_int && lower_limit_double) {
    std::cerr << "error: upper_limit and lower_limit must be the same type";
    throw;
  }

  std::vector<lapack_int> isuppz(2*dim+1);
  timer.start(timer_id::diagonalize_diagonalize);
  int info;
  if(mat.is_col_major())
    info = LAPACKE_dsyevr(LAPACK_COL_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, NULL, ldim_mat, &isuppz[0]);
  else
    info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, NULL, ldim_mat, &isuppz[0]);
  timer.stop(timer_id::diagonalize_diagonalize);
  timer.start(timer_id::diagonalize_finalize);
  if (info) {
    std::cerr << "error at dsyevr function. info=" << info << std::endl;
    exit(1);
  }
  params_out.set("m", m);
  params_out.set("isuppz", isuppz);
  
  if (params.get_bool("verbose")) {
    if (range == 'A')
      std::cout << "All eigenvalues/eigenvectors were requested" << std::endl;
    else if (upper_limit_double && lower_limit_double)
      std::cout << "Eigenvalues/eigenvectors contained in the interval [" << vl << ", " << vu << "]" << " were requested" << std::endl;
    else if (upper_limit_int && lower_limit_int)
      std::cout << "Eigenvalues/eigenvectors from " << il << "th" << " to " << iu << "th" << " were requested" << std::endl;
    std::cout << "The number of found eigenvalues are " << m << std::endl;
    std::cout << "The " << matrix_part << " part of the matrix was used" << std::endl;
    std::cout << "abstol=" << abstol << std::endl;
  }
  timer.stop(timer_id::diagonalize_finalize);
  return info;
}


// dsyevr eigenvalues / eigenvectors
template<typename MATRIX_MAJOR>
int diagonalize_r(localized_matrix<double, MATRIX_MAJOR>& mat, double* eigvals,
		  localized_matrix<double, MATRIX_MAJOR>& eigvecs,
		  rokko::parameters const& params, timer& timer) {
  char jobz = 'V';  // eigenvalues / eigenvectors

  rokko::parameters params_out;
  int dim = mat.outerSize();
  int ldim_mat = mat.innerSize();
  int ldim_eigvec = eigvecs.innerSize();
  std::vector<lapack_int> isuppz(2*dim+1);

  lapack_int m;  // output: found eigenvalues
  double abstol = 0.;  // defalut value = 0
  get_key(params, "abstol", abstol);

  std::string matrix_part = "upper"; // default is "upper"
  char uplow = 'U';
  if (params.defined("uplow"))
    matrix_part = params.get_string("uplow");
  if (params.defined("matrix_part"))
    matrix_part = params.get_string("matrix_part");
  if ((matrix_part == "upper") || (matrix_part == "U"))
    matrix_part = "upper";  uplow = 'U';
  if ((matrix_part == "lower") || (matrix_part == "L"))
    matrix_part = "lower";  uplow = 'L';

  char range = 'A';  // default is 'A'
  lapack_int il = 0, iu = 0;
  double vl = 0, vu = 0;
  bool upper_limit_double = get_key(params, "upper_limit", vu);
  bool upper_limit_int = get_key(params, "upper_limit", iu);
  bool lower_limit_double = get_key(params, "lower_limit", vl);
  bool lower_limit_int = get_key(params, "lower_limit", il);
  if (upper_limit_int && lower_limit_int)   range = 'I';
  if (upper_limit_double && lower_limit_double)   range = 'V';
  if (upper_limit_int && lower_limit_double) {
    std::cerr << "error: upper_limit and lower_limit must be the same type";
    throw;
  }

  timer.start(timer_id::diagonalize_diagonalize);
  int info;
  if(mat.is_col_major())
    info = LAPACKE_dsyevr(LAPACK_COL_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, &eigvecs(0,0), ldim_eigvec, &isuppz[0]);
  else
    info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, jobz, range, uplow, dim, &mat(0,0), ldim_mat, vl, vu, il, iu, abstol, &m, eigvals, &eigvecs(0,0), ldim_eigvec, &isuppz[0]);
  timer.stop(timer_id::diagonalize_diagonalize);
  timer.start(timer_id::diagonalize_finalize);
  if (info) {
    std::cerr << "error at dsyevr function. info=" << info << std::endl;
    exit(1);
  }
  params_out.set("m", m);
  params_out.set("isuppz", isuppz);
  
  if (params.get_bool("verbose")) {
    if (range == 'A')
      std::cout << "All eigenvalues/eigenvectors were requested" << std::endl;
    else if (upper_limit_double && lower_limit_double)
      std::cout << "Eigenvalues/eigenvectors contained in the interval [" << vl << ", " << vu << "]" << " were requested" << std::endl;
    else if (upper_limit_int && lower_limit_int)
      std::cout << "Eigenvalues/eigenvectors from " << il << "th" << " to " << iu << "th" << " were requested" << std::endl;
    std::cout << "The number of found eigenvalues are " << m << std::endl;
    std::cout << "The " << matrix_part << " part of the matrix was used" << std::endl;
    std::cout << "abstol=" << abstol << std::endl;
  }
  timer.stop(timer_id::diagonalize_finalize);
  return info;
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_R_HPP
