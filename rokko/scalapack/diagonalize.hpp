#ifndef ROKKO_SCALAPACK_DIAGONALIZE_H
#define ROKKO_SCALAPACK_DIAGONALIZE_H

#include <mpi.h>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/scalapack/blacs.hpp>
#include <rokko/scalapack/scalapack.hpp>

namespace rokko {
namespace scalapack {

template<typename MATRIX_MAJOR>
int diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, double* eigvals, distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  int ictxt;
  int info;

  const int ZERO=0, ONE=1;
  long MINUS_ONE = -1;
  blacs_get_(MINUS_ONE, ZERO, ictxt);

  char char_grid_major;
  if(mat.g.is_row_major())  char_grid_major = 'R';
  else  char_grid_major = 'C';

  blacs_gridinit_(ictxt, &char_grid_major, mat.g.get_nprow(), mat.g.get_npcol()); // ColがMPI_Comm_createと互換

  int dim = mat.get_m_global();
  int desc[9];
  descinit_(desc, mat.get_m_global(), mat.get_n_global(), mat.get_mb(), mat.get_nb(), ZERO, ZERO, ictxt, mat.get_lld(), info);
  if (info) {
    std::cerr << "error " << info << " at descinit function of descA " << "mA=" << mat.get_m_local() << "  nA=" << mat.get_n_local() << "  lld=" << mat.get_lld() << "." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 89);
  }

#ifdef DEBUG
  for (int proc=0; proc<mat.nprocs; ++proc) {
    if (proc == mat.g.get_myrank()) {
      std::cout << "pdsyev:proc=" << proc << " m_global=" << mat.get_m_global() << "  n_global=" << mat.get_n_global() << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  double* work = new double[1];
  long lwork = -1;

  // work配列のサイズの問い合わせ
  char* V = const_cast<char*>("V");
  char* U = const_cast<char*>("U");
  pdsyev_(V,  U,  dim,  mat.get_array_pointer(), ONE,  ONE,  desc, eigvals, eigvecs.get_array_pointer(), ONE, ONE,
 	   desc, work, lwork, info);

  lwork = work[0];
  delete[] work;
  work = new double [lwork];
  if (work == NULL) {
    std::cerr << "failed to allocate work. info=" << info << std::endl;
    return info;
  }

  // 固有値分解
  pdsyev_(V,  U,  dim,  mat.get_array_pointer(),  ONE,  ONE,  desc, eigvals, eigvecs.get_array_pointer(), ONE, ONE,
          desc, work, lwork, info);

  if (info) {
    std::cerr << "error at pdsyev function. info=" << info  << std::endl;
    exit(1);
  }

  delete[] work;
  return info;
}

template<class MATRIX_MAJOR>
int diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  return diagonalize(mat, &eigvals[0], eigvecs, timer_in);
}

} // namespace scalapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_DIAGONALIZE_H
