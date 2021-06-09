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
#include <rokko/scalapack/psyevx.hpp>
#include <rokko/utility/timer.hpp>

#include <mpi.h>

namespace rokko {
namespace scalapack {

// eigenvalues & eigenvectors
template<typename T, typename VEC>
parameters diagonalize_psyevx(distributed_matrix<T, rokko::matrix_col_major>& mat,
			       VEC& eigvals, distributed_matrix<T, rokko::matrix_col_major>& eigvecs,
			       parameters const& params) {
  parameters params_out;
  const char uplow = lapack::get_matrix_part(params);
  real_t<T> vl, vu;
  int il, iu;
  const char range = lapack::get_eigenvalues_range(params, vl, vu, il, iu);
  const int ictxt = mat.get_grid().get_blacs_context();
  real_t<T> abstol = params.defined("abstol") ? params.get<real_t<T>>("abstol") : cscalapack_pdlamch(ictxt, 'U');
  real_t<T> orfac = params.defined("orfac") ? params.get<real_t<T>>("orfac") : -1.;  // default value is 10^{-3} for a minus value.
  std::vector<int> ifail(mat.get_m_global());
  std::vector<int> iclustr(2 * mat.get_nprow() * mat.get_npcol());
  std::vector<real_t<T>> gap(mat.get_nprow() * mat.get_npcol());
  int m, nz;
  int info = psyevx(range, uplow, mat,
                    vl, vu, il, iu,
                    abstol, m, nz,
                    eigvals, orfac, eigvecs,
                    ifail, iclustr, gap);
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("nz", nz);
  params_out.set("ifail", ifail);
  params_out.set("iclustr", iclustr);
  params_out.set("gap", gap);
  if (params.get_bool("verbose")) {
    lapack::print_verbose("syevx", 'V', range, uplow, vl, vu, il, iu, params_out);
  }
  return params_out;
}

template<typename T, typename VEC>
parameters diagonalize_psyevx(distributed_matrix<T, rokko::matrix_row_major>& mat,
			       VEC& eigvals, distributed_matrix<T, rokko::matrix_row_major>& eigvecs,
			       parameters const& params) {
  throw std::invalid_argument("scalapack::diagonalize_psyevx() : scalapack doesn't support matrix_row_major.  Use scalapack with matrix_col_major.");
}

// only eigenvalues
template<typename T, typename VEC>
parameters diagonalize_psyevx(distributed_matrix<T, rokko::matrix_col_major>& mat,
			       VEC& eigvals,
			       parameters const& params) {
  rokko::parameters params_out;
  const char uplow = lapack::get_matrix_part(params);
  real_t<T> vl, vu;
  int il, iu;
  const char range = lapack::get_eigenvalues_range(params, vl, vu, il, iu);
  const int ictxt = mat.get_grid().get_blacs_context();
  real_t<T> abstol = params.defined("abstol") ? params.get<real_t<T>>("abstol") : cscalapack_pdlamch(ictxt, 'U');
  real_t<T> orfac = params.defined("orfac") ? params.get<real_t<T>>("orfac") : -1.;  // default value is 10^{-3} for a minus value.
  std::vector<int> ifail(mat.get_m_global());
  std::vector<int> iclustr(2 * mat.get_nprow() * mat.get_npcol());
  std::vector<real_t<T>> gap(mat.get_nprow() * mat.get_npcol());
  int m, nz;
  int info = psyevx(range, uplow, mat,
                    vl, vu, il, iu,
                    abstol, m, nz, eigvals, orfac,
                    ifail, iclustr, gap);
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("nz", nz);
  params_out.set("ifail", ifail);
  params_out.set("iclustr", iclustr);
  params_out.set("gap", gap);
  if (params.get_bool("verbose")) {
    lapack::print_verbose("syevx", 'N', range, uplow, vl, vu, il, iu, params_out);
  }
  return params_out;
}

template<typename T, typename VEC>
parameters diagonalize_psyevx(distributed_matrix<T, rokko::matrix_row_major>& mat,
			       VEC& eigvals,
			       parameters const& params) {
  throw std::invalid_argument("scalapack::diagonalize_psyevx() : scalapack doesn't support matrix_row_major.  Use scalapack with matrix_col_major.");
}

} // namespace scalapack
} // namespace rokko
