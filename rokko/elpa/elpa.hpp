#ifndef ROKKO_ELPA_HPP
#define ROKKO_ELPA_HPP

#include "mpi.h"

extern "C" {

void get_elpa_row_col_comms_wrap_(MPI_Fint* comm, int* my_prow, int* my_pcol, MPI_Fint* comm_rows,
                                  MPI_Fint* comm_cols);

void solve_evp_real_wrap_(int* na, int* nev, double* a, int* na_rows, double* ev, double* z,
                          int* nz_rows, int* nblk, MPI_Fint* comm_rows, MPI_Fint* comm_cols);

}

#endif // ROKKO_ELPA_HPP
