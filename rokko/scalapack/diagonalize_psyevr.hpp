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

#ifndef ROKKO_SCALAPACK_DIAGONALIZE_PSYEVR_HPP
#define ROKKO_SCALAPACK_DIAGONALIZE_PSYEVR_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/parameters.hpp>
#include <rokko/cscalapack.h>
#include <rokko/lapack/diagonalize_get_parameters.hpp>
#include <rokko/scalapack/psyevr.hpp>
#include <rokko/utility/timer.hpp>

#include <mpi.h>

namespace rokko {

namespace scalapack {

// eigenvalues & eigenvectors
template<typename T, typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_psyevr(distributed_matrix<T, MATRIX_MAJOR>& mat,
			VEC& eigvals, distributed_matrix<T, MATRIX_MAJOR>& eigvecs,
			parameters const& params) {
  parameters params_out;
  const char uplow = lapack::get_matrix_part(params);
  T vl = 0, vu = 0;
  int il = 0, iu = 0;
  const char range = lapack::get_eigenvalues_range(params, vu, vl, iu, il);

  int m, nz;
  int info = psyevr(range, uplow, mat,
                    vl, vu, il, iu, m, nz,
                    eigvals, eigvecs);
  if (info) {
    std::cerr << "error at pdsyevr function. info=" << info << std::endl;
  }
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("nz", nz);
  if ((mat.get_myrank() == 0) && params.get_bool("verbose")) {
    lapack::print_verbose("pdsyevr", 'V', range, uplow, vl, vu, il, iu, params_out);
  }
  return params_out;
}

// only eigenvalues
template<typename T, typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_psyevr(distributed_matrix<T, MATRIX_MAJOR>& mat,
			       VEC& eigvals,
			       parameters const& params) {
  parameters params_out;
  const char uplow = lapack::get_matrix_part(params);
  T vl = 0, vu = 0;
  int il = 0, iu = 0;
  const char range = lapack::get_eigenvalues_range(params, vu, vl, iu, il);
  int m, nz;
  int info = psyevr(range, uplow, mat,
                    vl, vu, il, iu, m, nz,
                    eigvals);

  if (info) {
    std::cerr << "error at pdsyevr function. info=" << info << std::endl;
  }
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("nz", nz);
  if ((mat.get_myrank() == 0) && params.get_bool("verbose")) {
    lapack::print_verbose("pdsyevr", 'N', range, uplow, vl, vu, il, iu, params_out);
  }
  return params_out;
}

} // namespace scalapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_DIAGONALIZE_PSYEVR_HPP