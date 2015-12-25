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
#include <rokko/elpa/diagonalize_get_parameters.hpp>
#include <rokko/utility/timer.hpp>

namespace rokko {
namespace elpa {

template<typename MATRIX_MAJOR>
parameters diagonalize_elpa2(distributed_matrix<double, MATRIX_MAJOR>& mat,
			     localized_vector<double>& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
			     parameters const& params) {
  parameters params_out;
  MPI_Fint comm_f = MPI_Comm_c2f(mat.get_grid().get_comm());
  int mpi_comm_rows, mpi_comm_cols;
  elpa_get_communicators(comm_f, mat.get_grid().get_myrow(), mat.get_grid().get_mycol(), &mpi_comm_rows, &mpi_comm_cols);
  
  // call eigenvalue routine
  int dim = mat.get_m_global();
  int nev = dim;
  get_nev(params, nev);
  //std::cout << "mat_lld=" << mat_lld << " mat_mb=" << mat_mb << std::endl;
  int kernel = ELPA2_REAL_KERNEL_GENERIC;
  get_kernel(params, kernel);
  int use_qr = 0;
  get_blocked_qr(params, use_qr);

  int info = elpa_solve_evp_real_2stage(dim, nev, mat.get_array_pointer(), mat.get_lld(), &eigvals[0],
					eigvecs.get_array_pointer(), eigvecs.get_lld(), mat.get_mb(), dim,
					mpi_comm_rows, mpi_comm_cols, comm_f, kernel, use_qr);
  params_out.set("info", info);
  return params_out;
}

template<typename MATRIX_MAJOR>
parameters diagonalize_elpa2(distributed_matrix<double, MATRIX_MAJOR>& mat,
			     localized_vector<double>& eigvals,
			     parameters const& params) {
  parameters params_out;
  std::cerr << "not yet implemented" << std::endl;
  throw;
  return params_out;
}

} // namespace elpa
} // namespace rokko

#endif // ROKKO_ELPA_DIAGONALIZE_ELPA2_HPP
