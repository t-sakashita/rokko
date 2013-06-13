#ifndef ROKKO_ELPA_HPP
#define ROKKO_ELPA_HPP

#include "mpi.h"

extern "C" {
  void get_elpa_row_col_comms_wrap_(const MPI_Comm& mpi_comm_world, const int& my_prow, const int& my_pcol, MPI_Fint& mpi_comm_rows, MPI_Fint& mpi_comm_cols);

  void solve_evp_real_wrap_(const int& na, const int& nev, double* a, const int& na_rows, double* ev, double* z, const int& nz_rows, const int& nblk, const MPI_Fint& mpi_comm_rows, const MPI_Fint& mpi_comm_cols);
}

#endif // ROKKO_ELPA_HPP
