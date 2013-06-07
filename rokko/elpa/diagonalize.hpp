#ifndef ROKKO_ELPA_DIAGONALIZE_H
#define ROKKO_ELPA_DIAGONALIZE_H

#include <mpi.h>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/elpa/blacs.hpp>
#include <rokko/elpa/elpa.hpp>

namespace rokko {
namespace elpa {

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

  MPI_Fint comm_f = (MPI_Fint) MPI_COMM_WORLD;
  MPI_Fint mpi_comm_rows_f, mpi_comm_cols_f;

  int nprow = mat.g.get_nprow();
  int npcol = mat.g.get_npcol();
 
  get_elpa_row_col_comms_wrap_(&comm_f, nprow, npcol,
			  &mpi_comm_rows_f, &mpi_comm_cols_f);

  int nblk = mat.get_nb();

  // 固有値分解
  timer_in.start(1);
  solve_evp_real_wrap_(dim, dim, mat.get_array_pointer(), dim, eigvals, eigvecs.get_array_pointer(), dim, nblk,
		 mpi_comm_rows_f, mpi_comm_cols_f);
  timer_in.stop(1);

  if (info) {
    std::cerr << "error at evp_real function. info=" << info  << std::endl;
    exit(1);
  }

  return info;
}

template<class MATRIX_MAJOR>
int diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  return diagonalize(mat, &eigvals[0], eigvecs, timer_in);
}

} // namespace elpa
} // namespace rokko

#endif // ROKKO_ELPA_DIAGONALIZE_H
