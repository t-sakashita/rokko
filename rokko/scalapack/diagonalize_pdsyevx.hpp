#ifndef ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVX_HPP
#define ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVX_HPP

#include <mpi.h>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/scalapack/blacs.hpp>
#include <rokko/scalapack/scalapack.hpp>

namespace rokko {
namespace scalapack {

template<typename MATRIX_MAJOR>
int diagonalize_x(distributed_matrix<MATRIX_MAJOR>& mat, double* eigvals, distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  int ictxt;
  int info;

  const int ZERO=0, ONE=1;
  long MINUS_ONE = -1;
  //blacs_pinfo_(mat.myrank, mat.nprocs);
  blacs_get_(MINUS_ONE, ZERO, ictxt);

  char char_grid_major;
  if(mat.get_grid().is_row_major())  char_grid_major = 'R';
  else  char_grid_major = 'C';

  //int tmp_nprow, tmp_npcol, tmp_myrow, tmp_mycol;  // we don't use these values
  blacs_gridinit_(ictxt, &char_grid_major, mat.get_grid().get_nprow(), mat.get_grid().get_npcol()); // ColがMPI_Comm_createと互換
  //blacs_gridinfo_(ictxt, tmp_nprow, tmp_npcol, tmp_myrow, tmp_mycol);

  int dim = mat.get_m_global();
  int desc[9];
  descinit_(desc, mat.get_m_global(), mat.get_n_global(), mat.get_mb(), mat.get_nb(), ZERO, ZERO, ictxt, mat.get_lld(), info);
  if (info) {
    std::cerr << "error " << info << " at descinit function of descA " << "mA=" << mat.get_m_local() << "  nA=" << mat.get_n_local() << "." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 89);
  }

#ifndef NDEBUG
  for (int proc=0; proc<mat.get_nprocs(); ++proc) {
    if (proc == mat.get_myrank()) {
      std::cout << "pdsyev:proc=" << proc << " m_global=" << mat.get_m_global() << "  n_global=" << mat.get_n_global() << "  mb=" << mat.get_mb() << "  nb=" << mat.get_nb() << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  int vl = 0, vu = 0;
  int il, iu;
  double abstol = pdlamch_(ictxt, 'U');  // get optimized absolute tolerance
  //std::cout << "abstol=" << abstol << std::endl;
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
  char* V = const_cast<char*>("V");
  char* A = const_cast<char*>("A");
  char* U = const_cast<char*>("U");

  timer_in.start(1);
  pdsyevx_(V, A, U, dim, mat.get_array_pointer(), ONE, ONE, desc,
           vl, vu, il, iu,
           abstol, num_eigval_found, num_eigvec_found, eigvals, orfac,
           eigvecs.get_array_pointer(), ONE, ONE, desc,
           work, lwork, iwork, liwork,

           ifail, iclustr, gap, info);
  timer_in.stop(1);

  if (info) {
    std::cerr << "error at pdsyevx function (query for sizes for workarrays." << std::endl;
    exit(1);
  }

  lwork = work[0];
  liwork = iwork[0];
  delete[] work;
  delete[] iwork;

  work = new double[lwork];
  iwork = new int[liwork];
  if (work == 0) {
    std::cerr << "failed to allocate work. info=" << info << std::endl;
    return info;
  }

  // 固有値分解
  pdsyevx_(V, A, U, dim, mat.get_array_pointer(), ONE, ONE, desc,
           vl, vu, il, iu,
           abstol, num_eigval_found, num_eigvec_found, eigvals, orfac,
           eigvecs.get_array_pointer(), ONE, ONE, desc,
           work, lwork, iwork, liwork,
           ifail, iclustr, gap, info);

  if (info) {
    std::cerr << "error at pdsyevx function. info=" << info << std::endl;
    exit(1);
  }

  delete[] work;
  delete[] iwork;
  delete[] ifail;
  delete[] iclustr;
  delete[] gap;
  return info;
}

template<class MATRIX_MAJOR>
int diagonalize_x(distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  return diagonalize_x(mat, &eigvals[0], eigvecs, timer_in);
}

} // namespace scalapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_DIAGONALIZE_PDSYEVX_HPP
