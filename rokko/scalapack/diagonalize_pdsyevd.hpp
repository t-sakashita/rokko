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

#ifndef ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVD_HPP
#define ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVD_HPP

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

// pdsyevd eigenvalues / eigenvectors
template<typename MATRIX_MAJOR>
int diagonalize_pdsyevd(distributed_matrix<double, MATRIX_MAJOR>& mat,
			localized_vector<double>& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			parameters const& params, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  char jobz = 'V';  // eigenvalues / eigenvectors
  std::string matrix_part = "upper"; // default is "upper"
  char uplow = 'U';
  lapack::get_matrix_part(params, matrix_part, uplow);

  int ictxt = ROKKO_blacs_get(-1, 0);
  char char_grid_major = blacs::set_grid_blacs(ictxt, mat);
  int dim = mat.get_m_global();
  int desc[9];
  blacs::set_desc(ictxt, mat, desc);
  int info;
  timer.stop(timer_id::diagonalize_initialize);

  timer.start(timer_id::diagonalize_diagonalize);
  info = ROKKO_pdsyevd(jobz, uplow, dim, mat.get_array_pointer(), 1, 1, desc, &eigvals[0],
		       eigvecs.get_array_pointer(), 1, 1, desc);
  timer.stop(timer_id::diagonalize_diagonalize);

  ROKKO_blacs_gridexit(&ictxt);
  return info;
}

} // namespace scalapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVD_HPP
