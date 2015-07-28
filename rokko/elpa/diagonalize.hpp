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

#ifndef ROKKO_ELPA_DIAGONALIZE_HPP
#define ROKKO_ELPA_DIAGONALIZE_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/parameters.hpp>
#include <rokko/elpa/elpa.hpp>
#include <rokko/utility/timer.hpp>

namespace rokko {
namespace elpa {

template<typename MATRIX_MAJOR>
parameters diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
		       localized_vector<double>& eigvals, distributed_matrix<double, MATRIX_MAJOR>& eigvecs,
		       parameters const& params) {
  parameters params_out;
  int dim = mat.get_m_global();
  MPI_Fint comm_f = MPI_Comm_c2f(mat.get_grid().get_comm());
  MPI_Fint comm_rows_f, comm_cols_f;
  int myrow = mat.get_grid().get_myrow();
  int mycol = mat.get_grid().get_mycol();
  get_elpa_row_col_comms_wrap_(&comm_f, &myrow, &mycol, &comm_rows_f, &comm_cols_f);

  // call eigenvalue routine
  int mat_lld = mat.get_lld();
  int mat_mb = mat.get_mb();
  int eigvecs_lld = eigvecs.get_lld();
  solve_evp_real_wrap_(&dim, &dim, mat.get_array_pointer(), &mat_lld, &eigvals[0],
		       eigvecs.get_array_pointer(), &eigvecs_lld, &mat_mb, &comm_rows_f, &comm_cols_f);
}

template<typename MATRIX_MAJOR>
parameters diagonalize(distributed_matrix<double, MATRIX_MAJOR>& mat,
		       localized_vector<double>& eigvals,
		       parameters const& params) {
  std::cerr << "not yet implemented" << std::endl;
  throw;
}

} // namespace elpa
} // namespace rokko

#endif // ROKKO_ELPA_DIAGONALIZE_HPP
