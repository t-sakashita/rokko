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

#ifndef ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVX_HPP
#define ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVX_HPP

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
template<typename T, typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_pdsyevx(distributed_matrix<T, MATRIX_MAJOR>& mat,
			       VEC& eigvals, distributed_matrix<T, MATRIX_MAJOR>& eigvecs,
			       parameters const& params) {
  parameters params_out;
  const char uplow = lapack::get_matrix_part(params);
  T vl, vu;
  int il, iu;
  const char range = lapack::get_eigenvalues_range(params, vl, vu, il, iu);
  const int ictxt = mat.get_grid().get_blacs_context();
  T abstol = params.defined("abstol") ? params.get<T>("abstol") : cscalapack_pdlamch(ictxt, 'U');
  T orfac = params.defined("orfac") ? params.get<T>("orfac") : -1.;  // default value is 10^{-3} for a minus value.
  std::vector<int> ifail(mat.get_m_global());
  std::vector<int> iclustr(2 * mat.get_nprow() * mat.get_npcol());
  std::vector<T> gap(mat.get_nprow() * mat.get_npcol());
  int m, nz;
  int info = psyevx(range, uplow, mat,
                    vl, vu, il, iu,
                    abstol, m, nz,
                    eigvals, orfac, eigvecs,
                    ifail.data(), iclustr.data(), gap.data());
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("nz", nz);
  params_out.set("ifail", ifail);
  params_out.set("iclustr", iclustr);
  params_out.set("gap", gap);
  if (params.get_bool("verbose")) {
    lapack::print_verbose("pdsyevx", 'V', range, uplow, vl, vu, il, iu, params_out);
  }
  return params_out;
}

// only eigenvalues
template<typename T, typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_pdsyevx(distributed_matrix<T, MATRIX_MAJOR>& mat,
			       VEC& eigvals,
			       parameters const& params) {
  rokko::parameters params_out;
  const char uplow = lapack::get_matrix_part(params);
  T vl, vu;
  int il, iu;
  const char range = lapack::get_eigenvalues_range(params, vl, vu, il, iu);
  const int ictxt = mat.get_grid().get_blacs_context();
  T abstol = params.defined("abstol") ? params.get<T>("abstol") : cscalapack_pdlamch(ictxt, 'U');
  T orfac = params.defined("orfac") ? params.get<T>("orfac") : -1.;  // default value is 10^{-3} for a minus value.
  std::vector<int> ifail(mat.get_m_global());
  std::vector<int> iclustr(2 * mat.get_nprow() * mat.get_npcol());
  std::vector<T> gap(mat.get_nprow() * mat.get_npcol());
  int m, nz;
  int info = psyevx(range, uplow, mat,
                    vl, vu, il, iu,
                    abstol, m, nz, eigvals, orfac,
                    ifail.data(), iclustr.data(), gap.data());
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("nz", nz);
  params_out.set("ifail", ifail);
  params_out.set("iclustr", iclustr);
  params_out.set("gap", gap);
  if (params.get_bool("verbose")) {
    lapack::print_verbose("pdsyevx", 'N', range, uplow, vl, vu, il, iu, params_out);
  }
  return params_out;
}

} // namespace scalapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVX_HPP
