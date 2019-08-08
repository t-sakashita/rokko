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

#ifndef ROKKO_ELPA_DIAGONALIZE_ELPA1_HPP
#define ROKKO_ELPA_DIAGONALIZE_ELPA1_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/parameters.hpp>
#include <rokko/elpa/elpa.h>
#include <rokko/elpa/diagonalize_get_parameters.hpp>
#include <rokko/utility/timer.hpp>

namespace rokko {
namespace elpa {

template<typename MATRIX_MAJOR>
parameters diagonalize_elpa1(distributed_matrix<double, MATRIX_MAJOR>& mat,
			     localized_vector<double>& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			     parameters const& params) {
  parameters params_out;
  if(mat.is_row_major())
    throw std::invalid_argument("elpa::diagonalize_elpa1() : elpa doesn't support matrix_row_major.  Use elpa with matrix_col_major.");
  MPI_Comm comm = mat.get_grid().get_comm();

  // call eigenvalue routine
  int dim = mat.get_m_global();
  int nev = dim;
  get_nev(params, nev);

  elpa_t handle;
  int error;
  handle = elpa_allocate(&error);

  /* Set parameters */
  elpa_set(handle, "na", dim, &error);
  assert_elpa_ok(error);
  
  elpa_set(handle, "nev", nev, &error);
  assert_elpa_ok(error);
  
  int na_rows = mat.get_m_local();
  elpa_set(handle, "local_nrows", na_rows, &error);
  assert_elpa_ok(error);
  
  int na_cols = mat.get_n_local();
  elpa_set(handle, "local_ncols", na_cols, &error);
  assert_elpa_ok(error);
  
  int nblk = mat.get_mb();
  elpa_set(handle, "nblk", nblk, &error);
  assert_elpa_ok(error);
  
  elpa_set(handle, "mpi_comm_parent", MPI_Comm_c2f(comm), &error);
  assert_elpa_ok(error);
  
  int my_prow = mat.get_grid().get_myrow();
  int my_pcol = mat.get_grid().get_mycol();
  MPI_Comm mpi_comm_rows, mpi_comm_cols;
  MPI_Comm_split(comm, my_pcol, my_prow, &mpi_comm_rows);  // color = my_pcol, key = my_prow
  MPI_Comm_split(comm, my_prow, my_pcol, &mpi_comm_cols);  // color = my_prow, key = my_pcol
  elpa_set(handle, "mpi_comm_rows", MPI_Comm_c2f(mpi_comm_rows), &error);
  assert_elpa_ok(error);
  elpa_set(handle, "mpi_comm_cols", MPI_Comm_c2f(mpi_comm_cols), &error);
  assert_elpa_ok(error);

  /* Setup */
  assert_elpa_ok(elpa_setup(handle));

  /* Set tunables */
  const int solver_enum = ELPA_SOLVER_1STAGE;
  elpa_set(handle, "solver", solver_enum, &error);
  assert_elpa_ok(error);
  
  int info;
  elpa_eigenvectors(handle, mat.get_array_pointer(), &eigvals[0], eigvecs.get_array_pointer(), &info);
  elpa_deallocate(handle, &error);

  params_out.set("info", info);
  return params_out;
}

template<typename MATRIX_MAJOR>
parameters diagonalize_elpa1(distributed_matrix<double, MATRIX_MAJOR>& mat,
			     localized_vector<double>& eigvals,
			     parameters const& params) {
  parameters params_out;
  if(mat.is_row_major())
    throw std::invalid_argument("elpa::diagonalize_elpa1() : elpa doesn't support matrix_row_major.  Use elpa with matrix_col_major.");
  MPI_Comm comm = mat.get_grid().get_comm();

  // call eigenvalue routine
  int dim = mat.get_m_global();
  int nev = 0;

  double *eigvecs = new double[dim*dim]; // No calculation of eigenvectors, but need as work array.

  elpa_t handle;
  int error;
  handle = elpa_allocate(&error);

  /* Set parameters */
  elpa_set(handle, "na", dim, &error);
  assert_elpa_ok(error);
  
  elpa_set(handle, "nev", nev, &error);
  assert_elpa_ok(error);
  
  int na_rows = mat.get_m_local();
  elpa_set(handle, "local_nrows", na_rows, &error);
  assert_elpa_ok(error);
  
  int na_cols = mat.get_n_local();
  elpa_set(handle, "local_ncols", na_cols, &error);
  assert_elpa_ok(error);
  
  int nblk = mat.get_mb();
  elpa_set(handle, "nblk", nblk, &error);
  assert_elpa_ok(error);
  
  elpa_set(handle, "mpi_comm_parent", MPI_Comm_c2f(comm), &error);
  assert_elpa_ok(error);
  
  int my_prow = mat.get_grid().get_myrow();
  int my_pcol = mat.get_grid().get_mycol();
  MPI_Comm mpi_comm_rows, mpi_comm_cols;
  MPI_Comm_split(comm, my_pcol, my_prow, &mpi_comm_rows);  // color = my_pcol, key = my_prow
  MPI_Comm_split(comm, my_prow, my_pcol, &mpi_comm_cols);  // color = my_prow, key = my_pcol
  elpa_set(handle, "mpi_comm_rows", MPI_Comm_c2f(mpi_comm_rows), &error);
  assert_elpa_ok(error);
  elpa_set(handle, "mpi_comm_cols", MPI_Comm_c2f(mpi_comm_cols), &error);
  assert_elpa_ok(error);
  
  /* Setup */
  assert_elpa_ok(elpa_setup(handle));

  /* Set tunables */
  const int solver_enum = ELPA_SOLVER_1STAGE;
  elpa_set(handle, "solver", solver_enum, &error);
  assert_elpa_ok(error);

  int info;
  elpa_eigenvalues(handle, mat.get_array_pointer(), &eigvals[0], &info);
  elpa_deallocate(handle, &error);

  delete[] eigvecs;
  params_out.set("info", info);
  return params_out;
}

} // namespace elpa
} // namespace rokko

#endif // ROKKO_ELPA_DIAGONALIZE_ELPA1_HPP
