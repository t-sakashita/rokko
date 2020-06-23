/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_ELPA_DIAGONALIZE_SET_PARAMETERS_HPP
#define ROKKO_ELPA_DIAGONALIZE_SET_PARAMETERS_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/parameters.hpp>
#include <rokko/elpa/elpa.h>
#include <rokko/elpa/diagonalize_get_parameters.hpp>

namespace rokko {
namespace elpa {

template<typename T, typename MATRIX_MAJOR>
void set_parameters(distributed_matrix<T, MATRIX_MAJOR> const& mat, parameters const& params, elpa_t handle) {
  int error;

  int nev = mat.get_m_global();
  get_nev(params, nev);

  elpa_set_integer(handle, "na", mat.get_m_global(), &error);
  assert_elpa_ok(error);

  elpa_set_integer(handle, "nev", nev, &error);
  assert_elpa_ok(error);

  int na_rows = mat.get_m_local();
  elpa_set_integer(handle, "local_nrows", na_rows, &error);
  assert_elpa_ok(error);

  int na_cols = mat.get_n_local();
  elpa_set_integer(handle, "local_ncols", na_cols, &error);
  assert_elpa_ok(error);

  int nblk = mat.get_mb();
  elpa_set_integer(handle, "nblk", nblk, &error);
  assert_elpa_ok(error);

  MPI_Comm comm = mat.get_grid().get_comm();
  elpa_set_integer(handle, "mpi_comm_parent", MPI_Comm_c2f(comm), &error);
  assert_elpa_ok(error);

  int my_prow = mat.get_grid().get_myrow();
  int my_pcol = mat.get_grid().get_mycol();
  MPI_Comm mpi_comm_rows, mpi_comm_cols;
  MPI_Comm_split(comm, my_pcol, my_prow, &mpi_comm_rows);  // color = my_pcol, key = my_prow
  MPI_Comm_split(comm, my_prow, my_pcol, &mpi_comm_cols);  // color = my_prow, key = my_pcol
  elpa_set_integer(handle, "mpi_comm_rows", MPI_Comm_c2f(mpi_comm_rows), &error);
  assert_elpa_ok(error);
  elpa_set_integer(handle, "mpi_comm_cols", MPI_Comm_c2f(mpi_comm_cols), &error);
  assert_elpa_ok(error);
}

} // namespace elpa
} // namespace rokko

#endif // ROKKO_ELPA_DIAGONALIZE_SET_PARAMETERS_HPP
