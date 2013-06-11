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
  if(mat.get_grid().is_row_major())  char_grid_major = 'R';
  else  char_grid_major = 'C';

  std::cout << "before: elpa_blacs" << std::endl;
  blacs_gridinit_(ictxt, &char_grid_major, mat.get_grid().get_nprow(), mat.get_grid().get_npcol()); // ColがMPI_Comm_createと互換
  std::cout << "after: elpa_blacs" << std::endl;

  int dim = mat.get_m_global();
  int desc[9];
  descinit_(desc, mat.get_m_global(), mat.get_n_global(), mat.get_mb(), mat.get_nb(), ZERO, ZERO, ictxt, mat.get_lld(), info);
  if (info) {
    std::cerr << "error " << info << " at descinit function of descA " << "mA=" << mat.get_m_local() << "  nA=" << mat.get_n_local() << "  lld=" << mat.get_lld() << "." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 89);
  }

  MPI_Fint comm_f = (MPI_Fint) MPI_COMM_WORLD;
  MPI_Fint mpi_comm_rows_f, mpi_comm_cols_f;

  int nprow = mat.get_grid().get_nprow();
  int npcol = mat.get_grid().get_npcol();

  std::cout << "elpa nprow=" << nprow << " npcol=" << npcol << std::endl;

  std::cout << "before: get_elpa_row_col_comms_wrap" << std::endl;

  get_elpa_row_col_comms_wrap_(&comm_f, mat.get_grid().get_myrow(), mat.get_grid().get_mycol(),
  			       &mpi_comm_rows_f, &mpi_comm_cols_f);
 
  int nblk = mat.get_mb();
  std::cout << "nblk=" << nblk << std::endl;

  double* mat_array = mat.get_array_pointer();
  double* eigvecs_array = eigvecs.get_array_pointer();

  std::cout << "before: solve_evp_real_wrap_" << std::endl;

  int m_local = mat.get_m_local();
  std::cout << "m_local=" << m_local << std::endl;

  // 固有値分解
  timer_in.start(1);
  //solve_evp_real_wrap2_(dim, dim, mat_array, m_local, eigvals, eigvecs_array, m_local, nblk, mat.get_grid().get_myrow(), mat.get_grid().get_mycol());    
  solve_evp_real_wrap_(dim, dim, mat_array, mat.get_m_local(), eigvals, eigvecs_array, eigvecs.get_m_local(), nblk,
  		       &mpi_comm_rows_f, &mpi_comm_cols_f);
  std::cout << "after: solve_evp_real_wrap_" << std::endl;
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
