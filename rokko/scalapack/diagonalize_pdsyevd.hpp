/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVD_HPP
#define ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVD_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/parameters.hpp>
#include <rokko/cscalapack.h>
#include <rokko/lapack/diagonalize_get_parameters.hpp>
#include <rokko/utility/timer.hpp>

#include <mpi.h>

namespace rokko {

namespace scalapack {

// eigenvalues & eigenvectors
template<typename T, typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_pdsyevd(distributed_matrix<T, MATRIX_MAJOR>& mat,
			       VEC& eigvals, distributed_matrix<T, MATRIX_MAJOR>& eigvecs,
			       parameters const& params) {
  parameters params_out;
  const char uplow = lapack::get_matrix_part(params);
  int info = psyevd(uplow, mat, eigvals, eigvecs);

  params_out.set("info", info);
  if (info) {
    std::cerr << "error at pdsyevd function. info=" << info << std::endl;
  }
  if ((mat.get_myrank() == 0) && params.get_bool("verbose")) {
    lapack::print_verbose("pdsyevd", 'V', uplow);
  }
  return params_out;
}

} // namespace scalapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVD_HPP
