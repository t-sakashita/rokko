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

#ifndef ROKKO_ELPA_DIAGONALIZE_ELPA2_HPP
#define ROKKO_ELPA_DIAGONALIZE_ELPA2_HPP

#include <rokko/distributed_matrix.hpp>
#include <rokko/parameters.hpp>
#include <rokko/elpa/elpa.h>
#include <rokko/elpa.hpp>
#include <rokko/elpa/diagonalize_get_parameters.hpp>
#include <rokko/elpa/diagonalize_set_parameters.hpp>
#include <rokko/utility/timer.hpp>

namespace rokko {
namespace elpa {

template<typename T, typename VEC>
parameters diagonalize_elpa2(distributed_matrix<T, rokko::matrix_col_major>& mat,
			     VEC& eigvals, distributed_matrix<T, rokko::matrix_col_major>& eigvecs,
			     parameters const& params) {
  parameters params_out;

  int error;
  elpa_t handle = elpa_allocate(&error);

  /* Set parameters */
  set_parameters(mat, params, handle);

  assert_elpa_ok( set(handle, "qr", get_blocked_qr(params)) );
  
  /* Setup */
  assert_elpa_ok(elpa_setup(handle));

  /* Set tunables */
  constexpr int solver_enum = ELPA_SOLVER_2STAGE;
  assert_elpa_ok( set(handle, "solver", solver_enum) );

  assert_elpa_ok( set(handle, "real_kernel", get_kernel(params)) );

  // call eigenvalue routine
  int info = elpa::diag(handle, mat, eigvals, eigvecs);
  assert_elpa_ok( deallocate(handle) );

  params_out.set("info", info);
  return params_out;
}

template<typename T, typename VEC>
parameters diagonalize_elpa2(distributed_matrix<T, rokko::matrix_row_major>& mat,
			     VEC& eigvals, distributed_matrix<T, rokko::matrix_row_major>& eigvecs,
			     parameters const& params) {
  throw std::invalid_argument("elpa::diagonalize_elpa2() : elpa doesn't support matrix_row_major.  Use elpa with matrix_col_major.");
}

template<typename T, typename VEC>
parameters diagonalize_elpa2(distributed_matrix<T, rokko::matrix_col_major>& mat,
			     VEC& eigvals,
			     parameters const& params) {
  parameters params_out;

  int error;
  elpa_t handle = elpa_allocate(&error);

  /* Set parameters */
  set_parameters(mat, params, handle);

  assert_elpa_ok( set(handle, "qr", get_blocked_qr(params)) );
  
  /* Setup */
  assert_elpa_ok(elpa_setup(handle));

  /* Set tunables */
  constexpr int solver_enum = ELPA_SOLVER_2STAGE;
  assert_elpa_ok( set(handle, "solver", solver_enum) );

  assert_elpa_ok( set(handle, "real_kernel", get_kernel(params)) );

  // call eigenvalue routine
  int info = elpa::diag(handle, mat, eigvals);
  assert_elpa_ok( deallocate(handle) );

  params_out.set("info", info);
  return params_out;
}

template<typename T, typename VEC>
parameters diagonalize_elpa2(distributed_matrix<T, rokko::matrix_row_major>& mat,
			     VEC& eigvals,
			     parameters const& params) {
  throw std::invalid_argument("elpa::diagonalize_elpa2() : elpa doesn't support matrix_row_major.  Use elpa with matrix_col_major.");
}

} // namespace elpa
} // namespace rokko

#endif // ROKKO_ELPA_DIAGONALIZE_ELPA2_HPP
