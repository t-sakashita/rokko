/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SCALAPACK_DIAGONALIZE_H
#define ROKKO_SCALAPACK_DIAGONALIZE_H

#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/blacs/blacs.h>
#include <rokko/scalapack/scalapack.h>
#include <rokko/utility/timer.hpp>

#include <mpi.h>

namespace rokko {
namespace scalapack {

template<typename MATRIX_MAJOR>
int diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
  distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  int ictxt;
  int info;

  ROKKO_blacs_get(-1, 0, &ictxt);

  char char_grid_major;
  if(mat.get_grid().is_row_major())  char_grid_major = 'R';
  else  char_grid_major = 'C';

  ROKKO_blacs_gridinit(&ictxt, char_grid_major, mat.get_grid().get_nprow(),
    mat.get_grid().get_npcol());
  int dim = mat.get_m_global();
  int desc[9];
  ROKKO_descinit(desc, mat.get_m_global(), mat.get_n_global(), mat.get_mb(), mat.get_nb(),
    0, 0, ictxt, mat.get_lld(), &info);
  if (info) {
    std::cerr << "error " << info << " at descinit function of descA " << "mA="
              << mat.get_m_local() << "  nA=" << mat.get_n_local() << "  lld=" << mat.get_lld()
              << "." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 89);
  }

  double* work = new double[1];
  long lwork = -1;

  // work配列のサイズの問い合わせ
  ROKKO_pdsyev('V', 'U', dim, mat.get_array_pointer(), 1, 1, desc, &eigvals[0],
    eigvecs.get_array_pointer(), 1, 1, desc, work, lwork, &info);

  lwork = work[0];
  delete[] work;
  work = new double [lwork];
  if (work == 0) {
    std::cerr << "failed to allocate work. info=" << info << std::endl;
    return info;
  }
  timer.stop(timer_id::diagonalize_initialize);

  // 固有値分解
  timer.start(timer_id::diagonalize_diagonalize);
  ROKKO_pdsyev('V', 'U', dim, mat.get_array_pointer(), 1, 1, desc, &eigvals[0],
    eigvecs.get_array_pointer(), 1, 1, desc, work, lwork, &info);
  if (info) {
    std::cerr << "error at pdsyev function. info=" << info  << std::endl;
    exit(1);
  }
  timer.stop(timer_id::diagonalize_diagonalize);

  timer.start(timer_id::diagonalize_finalize);
  delete[] work;
  ROKKO_blacs_gridexit(&ictxt);
  timer.stop(timer_id::diagonalize_finalize);
  return info;
}

} // namespace scalapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_DIAGONALIZE_H
