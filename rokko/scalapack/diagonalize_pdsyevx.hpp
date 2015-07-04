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

#ifndef ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVX_HPP
#define ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVX_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/parameters.hpp>
#include <rokko/blacs/blacs_wrap.h>
#include <rokko/blacs/utility_routines.hpp>
#include <rokko/scalapack/scalapack_wrap.h>
#include <rokko/lapack/diagonalize_get_parameters.hpp>
#include <rokko/utility/timer.hpp>

#include <mpi.h>

namespace rokko {
namespace scalapack {

template<typename MATRIX_MAJOR>
int diagonalize_pdsyevx(distributed_matrix<double, MATRIX_MAJOR>& mat,
			localized_vector<double>& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			parameters const& params, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  char jobz = 'V';  // eigenvalues / eigenvectors
  std::string matrix_part = "upper"; // default is "upper"
  char uplow = 'U';
  lapack::get_matrix_part(params, matrix_part, uplow);
  char range = 'A';  // default is 'A'
  double vl = 0, vu = 0;
  int il = 0, iu = 0;
  bool is_upper_value, is_upper_index, is_lower_value, is_lower_index;
  lapack::get_eigenvalues_range(params, matrix_part, range,
				vu, vl, iu, il,
				is_upper_value, is_lower_value, is_upper_index, is_lower_index);

  int ictxt = ROKKO_blacs_get(-1, 0);
  char char_grid_major = rokko::blacs::set_grid_blacs(ictxt, mat);
  int dim = mat.get_m_global();
  int desc[9];
  rokko::blacs::set_desc(ictxt, mat, desc);
  int m, nz;
  int info;
 
  double abstol = ROKKO_pdlamch(ictxt, 'U');
  //double abstol = 0.;  // defalut value = 0
  //get_key(params, "abstol", abstol);

  int num_eigval_found, num_eigvec_found;
  int orfac = -1;  // default value 10^{-3} is used.
  int* ifail = new int[dim];
  int* iclustr = new int[2 * mat.get_nprow() * mat.get_npcol()];
  double* gap = new double[mat.get_nprow() * mat.get_npcol()];
  if (ifail == 0) {
    std::cerr << "failed to allocate ifail." << std::endl;
    exit(1);
  }
  if (iclustr == 0) {
    std::cerr << "failed to allocate iclustr." << std::endl;
    exit(1);
  }
  if (gap == 0) {
    std::cerr << "failed to allocate gap." << std::endl;
    exit(1);
  }
  timer.stop(timer_id::diagonalize_initialize);

  timer.start(timer_id::diagonalize_diagonalize);
  info = ROKKO_pdsyevx(jobz, range, uplow, dim, mat.get_array_pointer(), 1, 1, desc, vl, vu, il, iu,
		       abstol, &num_eigval_found, &num_eigvec_found, &eigvals[0], orfac,
		       eigvecs.get_array_pointer(), 1, 1, desc,
		       ifail, iclustr, gap);
  timer.stop(timer_id::diagonalize_diagonalize);

  timer.start(timer_id::diagonalize_finalize);
  delete[] ifail;
  delete[] iclustr;
  delete[] gap;
  ROKKO_blacs_gridexit(&ictxt);
  timer.stop(timer_id::diagonalize_finalize);
  return info;
}

} // namespace scalapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVX_HPP
