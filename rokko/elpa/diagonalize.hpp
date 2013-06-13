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
  int dim = mat.get_m_global();

  MPI_Fint comm_f = MPI_COMM_WORLD;
  MPI_Fint mpi_comm_rows_f, mpi_comm_cols_f;

  get_elpa_row_col_comms_wrap_(comm_f, mat.get_grid().get_myrow(), mat.get_grid().get_mycol(),
  			       mpi_comm_rows_f, mpi_comm_cols_f);

  // call eigenvalue routine
  timer_in.start(1);
  solve_evp_real_wrap_(dim, dim, mat.get_array_pointer(), mat.get_lld(), eigvals, eigvecs.get_array_pointer(), eigvecs.get_lld(), mat.get_mb(),
  		       mpi_comm_rows_f, mpi_comm_cols_f);
  timer_in.stop(1);
}

template<class MATRIX_MAJOR>
int diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  return diagonalize(mat, &eigvals[0], eigvecs, timer_in);
}

} // namespace elpa
} // namespace rokko

#endif // ROKKO_ELPA_DIAGONALIZE_H
