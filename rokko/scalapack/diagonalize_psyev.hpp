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

#ifndef ROKKO_SCALAPACK_DIAGONALIZE_PSYEV_HPP
#define ROKKO_SCALAPACK_DIAGONALIZE_PSYEV_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/parameters.hpp>
#include <rokko/scalapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>
#include <rokko/utility/timer.hpp>

namespace rokko {
namespace scalapack {

// eigenvalues & eigenvectors
template<typename T, typename VEC>
parameters diagonalize_psyev(distributed_matrix<T, rokko::matrix_col_major>& mat,
			      VEC& eigvals, distributed_matrix<T, rokko::matrix_col_major>& eigvecs,
			      parameters const& params) {
  parameters params_out;
  const char uplow = lapack::get_matrix_part(params);
  int info = psyev(uplow, mat, eigvals, eigvecs);
  params_out.set("info", info);
  if (info) {
    std::cerr << "error at pdsyev function. info=" << info << std::endl;
  }
  if ((mat.get_myrank() == 0) && params.get_bool("verbose")) {
    lapack::print_verbose("syev", 'V', uplow);
  }
  return params_out;
}

template<typename T, typename VEC>
parameters diagonalize_psyev(distributed_matrix<T, rokko::matrix_row_major>& mat,
			      VEC& eigvals, distributed_matrix<T, rokko::matrix_row_major>& eigvecs,
			      parameters const& params) {
  throw std::invalid_argument("scalapack::diagonalize_psyev() : scalapack doesn't support matrix_row_major.  Use scalapack with matrix_col_major.");
}

// only eigenvalues
template<typename T, typename VEC>
parameters diagonalize_psyev(distributed_matrix<T, rokko::matrix_col_major>& mat,
			      VEC& eigvals,
			      parameters const& params) {
  parameters params_out;
  const char uplow = lapack::get_matrix_part(params);
  int info = psyev(uplow, mat, eigvals);
  params_out.set("info", info);
  if (info) {
    std::cerr << "error at pdsyev function. info=" << info << std::endl;
  }
  if ((mat.get_myrank() == 0) && params.get_bool("verbose")) {
    lapack::print_verbose("syev", 'N', uplow);
  }
  return params_out;
}

template<typename T, typename VEC>
parameters diagonalize_psyev(distributed_matrix<T, rokko::matrix_row_major>& mat,
			      VEC& eigvals,
			      parameters const& params) {
  throw std::invalid_argument("scalapack::diagonalize_psyev() : scalapack doesn't support matrix_row_major.  Use scalapack with matrix_col_major.");
}

} // namespace scalapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_DIAGONALIZE_PSYEV_HPP
