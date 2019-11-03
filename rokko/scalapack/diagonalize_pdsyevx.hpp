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
#include <rokko/utility/timer.hpp>

#include <mpi.h>

namespace rokko {
namespace scalapack {

// pdsyevx eigenvalues / eigenvectors
template<typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_pdsyevx(distributed_matrix<double, MATRIX_MAJOR>& mat,
			       VEC& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			       parameters const& params) {
  parameters params_out;
  const char jobz = 'V';  // eigenvalues / eigenvectors
  const char uplow = lapack::get_matrix_part(params);
  double vl, vu;
  int il, iu;
  const char range = lapack::get_eigenvalues_range(params, vl, vu, il, iu);
  const int ictxt = mat.get_grid().get_blacs_context();
  const int* desc = mat.get_mapping().get_blacs_descriptor();
  double abstol = cscalapack_pdlamch(ictxt, 'U');
  //get_key(params, "abstol", abstol);
  int orfac = -1;  // default value 10^{-3} is used.
  std::vector<int> ifail(mat.get_m_global());
  std::vector<int> iclustr(2 * mat.get_nprow() * mat.get_npcol());
  std::vector<double> gap(mat.get_nprow() * mat.get_npcol());
  int m, nz;
  int info = cscalapack_pdsyevx(jobz, range, uplow, mat.get_m_global(), mat.get_array_pointer(), 0, 0, desc, vl, vu, il, iu,
		       abstol, &m, &nz, &eigvals[0], orfac,
		       eigvecs.get_array_pointer(), 0, 0, desc,
		       &ifail[0], &iclustr[0], &gap[0]);
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("nz", nz);
  params_out.set("ifail", ifail);
  params_out.set("iclustr", iclustr);
  params_out.set("gap", gap);
  if (params.get_bool("verbose")) {
    lapack::print_verbose("pdsyevx", jobz, range, uplow, vl, vu, il, iu, params_out);
  }
  return params_out;
}

// pdsyevx only eigenvalues
template<typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_pdsyevx(distributed_matrix<double, MATRIX_MAJOR>& mat,
			       VEC& eigvals,
			       parameters const& params) {
  rokko::parameters params_out;
  const char jobz = 'N';  // only eigenvalues
  const char uplow = lapack::get_matrix_part(params);
  double vl, vu;
  int il, iu;
  const char range = lapack::get_eigenvalues_range(params, vl, vu, il, iu);
  const int ictxt = mat.get_grid().get_blacs_context();
  const int* desc = mat.get_mapping().get_blacs_descriptor();
  double abstol = cscalapack_pdlamch(ictxt, 'U');
  //get_key(params, "abstol", abstol);
  int orfac = -1;  // default value 10^{-3} is used.
  std::vector<int> ifail(mat.get_m_global());
  std::vector<int> iclustr(2 * mat.get_nprow() * mat.get_npcol());
  std::vector<double> gap(mat.get_nprow() * mat.get_npcol());
  int m, nz;
  int info = cscalapack_pdsyevx(jobz, range, uplow, mat.get_m_global(), mat.get_array_pointer(), 0, 0, desc, vl, vu, il, iu,
		       abstol, &m, &nz, &eigvals[0], orfac,
		       NULL, 0, 0, desc,
		       &ifail[0], &iclustr[0], &gap[0]);
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("nz", nz);
  params_out.set("ifail", ifail);
  params_out.set("iclustr", iclustr);
  params_out.set("gap", gap);
  if (params.get_bool("verbose")) {
    lapack::print_verbose("pdsyevx", jobz, range, uplow, vl, vu, il, iu, params_out);
  }
  return params_out;
}

} // namespace scalapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVX_HPP
