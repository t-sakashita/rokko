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
  blacs_pinfo_(mat.myrank, mat.nprocs);
  blacs_get_(MINUS_ONE, ZERO, ictxt);

  char char_grid_major;
  if(mat.g.is_row_major())  char_grid_major = 'R';
  else  char_grid_major = 'C';

  blacs_gridinit_(ictxt, &char_grid_major, mat.nprow, mat.npcol); // ColがMPI_Comm_createと互換
  blacs_gridinfo_(ictxt, mat.nprow, mat.npcol, mat.myrow, mat.mycol);

  int dim = mat.m_global;
  //std::cout << "pdsyev_dim=" << dim << std::endl;
  int desc[9];

  int lld = mat.m_local;
  if (lld == 0) lld = 1;
  descinit_(desc, mat.m_global, mat.n_global, mat.mb, mat.nb, ZERO, ZERO, ictxt, lld, info);
  if (info) {
    std::cerr << "error " << info << " at descinit function of descA " << "mA=" << mat.m_local << "  nA=" << mat.n_local << "  lld=" << lld << "." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 89);
  }

#ifdef DEBUG
  for (int proc=0; proc<mat.nprocs; ++proc) {
    if (proc == mat.myrank) {
      std::cout << "pdsyev:proc=" << proc << " m_global=" << mat.m_global << "  n_global=" << mat.n_global << "  mb=" << mat.mb << "  nb=" << mat.nb << " lld=" << lld << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  double* work = new double[1];
  long lwork = -1;

  // work配列のサイズの問い合わせ
  char* V = const_cast<char*>("V");
  char* U = const_cast<char*>("U");
  pdsyev_( V,  U,  dim,  mat.array, ONE,  ONE,  desc, eigvals, eigvecs.array, ONE, ONE,
 	   desc, work, lwork, info );

  lwork = work[0];
  delete[] work;
  work = new double [lwork];
  if (work == NULL) {
    std::cerr << "failed to allocate work. info=" << info << std::endl;
    return info;
  }

  // 固有値分解
  pdsyev_( V,  U,  dim,  mat.array,  ONE,  ONE,  desc, eigvals, eigvecs.array, ONE, ONE,
  	   desc, work, lwork, info );

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
