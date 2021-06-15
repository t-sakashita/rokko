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

#pragma once

#include <rokko/distributed_matrix.hpp>
#include <rokko/parameters.hpp>
#include <rokko/cscalapack.h>
#include <rokko/lapack/diagonalize_get_parameters.hpp>
#include <rokko/utility/timer.hpp>

#include <mpi.h>

namespace rokko {

namespace scalapack {

// eigenvalues & eigenvectors
template<typename T, typename VEC>
parameters diagonalize_psyevd(distributed_matrix<T, rokko::matrix_col_major>& mat,
			       VEC& eigvals, distributed_matrix<T, rokko::matrix_col_major>& eigvecs,
			       parameters const& params) {
  parameters params_out;
  const char uplow = lapack::get_matrix_part(params);
  int info = psyevd(uplow, mat, eigvals, eigvecs);

  params_out.set("info", info);
  if (info) {
    std::cerr << "error at pdsyevd function. info=" << info << std::endl;
  }
  if ((mat.get_myrank() == 0) && params.get_bool("verbose")) {
    lapack::print_verbose("syevd", 'V', uplow);
  }
  return params_out;
}

template<typename T, typename VEC>
parameters diagonalize_psyevd(distributed_matrix<T, rokko::matrix_row_major>& /* mat */,
                              VEC& /* eigvals */, distributed_matrix<T, rokko::matrix_row_major>& /* eigvecs */,
                              parameters const& /* params */) {
  throw std::invalid_argument("scalapack::diagonalize_psyevd() : scalapack doesn't support matrix_row_major.  Use scalapack with matrix_col_major.");
}

} // namespace scalapack
} // namespace rokko
