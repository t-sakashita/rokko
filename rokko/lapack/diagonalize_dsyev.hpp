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

#ifndef ROKKO_LAPACK_DIAGONALIZE_DSYEV_HPP
#define ROKKO_LAPACK_DIAGONALIZE_DSYEV_HPP

#include <rokko/parameters.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>

namespace rokko {
namespace lapack {

// dsyev only eigenvalues
template<typename MATRIX_MAJOR>
parameters diagonalize_dsyev(localized_matrix<double, MATRIX_MAJOR>& mat, double* eigvals,
			     parameters const& params) {
  rokko::parameters params_out;
  char jobz = 'N';  // only eigenvalues
  char uplow = get_matrix_part(params);

  int dim = mat.outerSize();
  int ldim = mat.innerSize();
  int info;

  if(mat.is_col_major())
    info = LAPACKE_dsyev(LAPACK_COL_MAJOR, jobz, uplow, dim, &mat(0,0), ldim, eigvals);
  else
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, uplow, dim, &mat(0,0), ldim, eigvals);

  params_out.set("info", info);
  if (info) {
    std::cerr << "error at dsyev function. info=" << info  << std::endl;
    exit(1);
  }
  if (params.get_bool("verbose")) {
    print_verbose("dsyev", jobz, uplow);
  }

  return params_out;
}

// dsyev eigenvalues / eigenvectors
template<typename MATRIX_MAJOR>
parameters diagonalize_dsyev(localized_matrix<double, MATRIX_MAJOR>& mat, double* eigvals,
			     localized_matrix<double, MATRIX_MAJOR>& eigvecs,
			     parameters const& params) {
  rokko::parameters params_out;
  char jobz = 'V';  // eigenvalues / eigenvectors
  char uplow = get_matrix_part(params);

  int dim = mat.outerSize();
  int ldim = mat.innerSize();
  int info;

  if(mat.is_col_major())
    info = LAPACKE_dsyev(LAPACK_COL_MAJOR, jobz, uplow, dim, &mat(0,0), ldim, eigvals);
  else
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, uplow, dim, &mat(0,0), ldim, eigvals);

  eigvecs = mat;
  params_out.set("info", info);
  if (info) {
    std::cerr << "error at dsyev function. info=" << info  << std::endl;
    exit(1);
  }
  if (params.get_bool("verbose")) {
    print_verbose("dsyev", jobz, uplow);
  }

  return params_out;
}

} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_DSYEV_HPP
