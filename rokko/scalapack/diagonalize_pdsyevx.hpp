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
#include <rokko/blacs/blacs_wrap.h>
#include <rokko/scalapack/scalapack_wrap.h>
#include <rokko/utility/timer.hpp>

#include <mpi.h>

namespace rokko {
namespace scalapack {

template<typename MATRIX_MAJOR>
int diagonalize_x(distributed_matrix<double, MATRIX_MAJOR>& mat, localized_vector<double>& eigvals,
  distributed_matrix<double, MATRIX_MAJOR>& eigvecs, timer& timer) {
  timer.start(timer_id::diagonalize_initialize);
  
  int ictxt = ROKKO_blacs_get(-1, 0);

  char char_grid_major;
  if(mat.get_grid().is_row_major())  char_grid_major = 'R';
  else  char_grid_major = 'C';

  ROKKO_blacs_gridinit(&ictxt, char_grid_major, mat.get_grid().get_nprow(),
    mat.get_grid().get_npcol());
  int dim = mat.get_m_global();
  int desc[9];
  int info = ROKKO_descinit(desc, mat.get_m_global(), mat.get_n_global(), mat.get_mb(),
                            mat.get_nb(), 0, 0, ictxt, mat.get_lld());
  if (info) {
    std::cerr << "error " << info << " at descinit function of descA " << "mA="
              << mat.get_m_local() << "  nA=" << mat.get_n_local() << "  lld=" << mat.get_lld()
              << "." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 89);
  }

  int vl = 0, vu = 0;
  int il, iu;
  double abstol = ROKKO_pdlamch(ictxt, 'U');
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

  double* work = new double[1];
  int* iwork = new int[1];
  long lwork = -1;
  long liwork = -1;

  // work配列のサイズの問い合わせ
  info = ROKKO_pdsyevx('V', 'A', 'U', dim, mat.get_array_pointer(), 1, 1, desc, vl, vu, il, iu,
    abstol, &num_eigval_found, &num_eigvec_found, &eigvals[0], orfac,
    eigvecs.get_array_pointer(), 1, 1, desc,
    work, lwork, iwork, liwork, ifail, iclustr, gap);

  lwork = work[0];
  delete[] work;
  work = new double[lwork];
  liwork = iwork[0];
  delete[] iwork;
  iwork = new int[liwork];
  if (work == 0 || iwork == 0) {
    std::cerr << "failed to allocate work. info=" << info << std::endl;
    return info;
  }
  timer.stop(timer_id::diagonalize_initialize);

  // 固有値分解
  timer.start(timer_id::diagonalize_diagonalize);
  info = ROKKO_pdsyevx('V', 'A', 'U', dim, mat.get_array_pointer(), 1, 1, desc, vl, vu, il, iu,
    abstol, &num_eigval_found, &num_eigvec_found, &eigvals[0], orfac,
    eigvecs.get_array_pointer(), 1, 1, desc,
    work, lwork, iwork, liwork, ifail, iclustr, gap);
  if (info) {
    std::cerr << "error at pdsyevx function. info=" << info << std::endl;
    exit(1);
  }
  timer.stop(timer_id::diagonalize_diagonalize);

  timer.start(timer_id::diagonalize_finalize);
  delete[] work;
  delete[] iwork;
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
