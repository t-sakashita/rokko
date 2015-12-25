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

#ifndef ROKKO_ELPA_DIAGONALIZE_ELPA2_HPP
#define ROKKO_ELPA_DIAGONALIZE_ELPA2_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/parameters.hpp>
#include <rokko/elpa/elpa.h>
#include <rokko/utility/timer.hpp>

namespace rokko {
namespace elpa {

template<typename MATRIX_MAJOR>
parameters diagonalize_elpa2(distributed_matrix<double, MATRIX_MAJOR>& mat,
			     localized_vector<double>& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			     parameters const& params) {
  parameters params_out;
  int dim = mat.get_m_global();
  MPI_Fint comm_f = MPI_Comm_c2f(mat.get_grid().get_comm());
  MPI_Fint comm_rows_f, comm_cols_f;
  int mpi_comm_rows, mpi_comm_cols;
  int myrow = mat.get_grid().get_myrow();
  int mycol = mat.get_grid().get_mycol();
  elpa_get_communicators(comm_f, myrow, mycol, &mpi_comm_rows, &mpi_comm_cols);
  
  // call eigenvalue routine
  int mat_lld = mat.get_lld();
  int mat_mb = mat.get_mb();
  int eigvecs_lld = eigvecs.get_lld();
  //std::cout << "mat_lld=" << mat_lld << " mat_mb=" << mat_mb << std::endl;
  int use_qr = 0;
  int THIS_REAL_ELPA_KERNEL_API = ELPA2_REAL_KERNEL_GENERIC;
  elpa_solve_evp_real_2stage(dim, dim, mat.get_array_pointer(), mat_lld, &eigvals[0],
  			     eigvecs.get_array_pointer(), eigvecs_lld, mat_mb, dim,
			     mpi_comm_rows, mpi_comm_cols, comm_f, THIS_REAL_ELPA_KERNEL_API, use_qr);

  return params;
}

template<typename MATRIX_MAJOR>
parameters diagonalize_elpa2(distributed_matrix<double, MATRIX_MAJOR>& mat,
			     localized_vector<double>& eigvals,
			     parameters const& params) {
  parameters params_out;
  std::cerr << "not yet implemented" << std::endl;
  throw;
  return params;
}

} // namespace elpa
} // namespace rokko

#endif // ROKKO_ELPA_DIAGONALIZE_ELPA2_HPP
