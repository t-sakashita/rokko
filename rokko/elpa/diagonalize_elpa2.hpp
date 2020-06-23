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

template<typename T, typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_elpa2(distributed_matrix<T, MATRIX_MAJOR>& mat,
			     VEC& eigvals, distributed_matrix<T, MATRIX_MAJOR>& eigvecs,
			     parameters const& params) {
  parameters params_out;
  if(mat.is_row_major())
    throw std::invalid_argument("elpa::diagonalize_elpa2() : elpa doesn't support matrix_row_major.  Use elpa with matrix_col_major.");

  elpa_t handle;
  int error;
  handle = elpa_allocate(&error);

  /* Set parameters */
  set_parameters(mat, params, handle);

  int use_qr = 0;
  get_blocked_qr(params, use_qr);
  elpa_set_integer(handle, "qr", use_qr, &error);
  
  /* Setup */
  assert_elpa_ok(elpa_setup(handle));

  /* Set tunables */
  constexpr int solver_enum = ELPA_SOLVER_2STAGE;
  elpa_set_integer(handle, "solver", solver_enum, &error);
  assert_elpa_ok(error);

  int kernel = get_kernel(params);
  elpa_set_integer(handle, "real_kernel", kernel, &error);
  assert_elpa_ok(error);

  // call eigenvalue routine
  int info = elpa::diag(handle, mat, eigvals, eigvecs);
  elpa_deallocate(handle, &error);

  params_out.set("info", info);
  return params_out;
}

template<typename T, typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_elpa2(distributed_matrix<T, MATRIX_MAJOR>& mat,
			     VEC& eigvals,
			     parameters const& params) {
  parameters params_out;
  
  elpa_t handle;
  int error;
  handle = elpa_allocate(&error);

  /* Set parameters */
  set_parameters(mat, params, handle);

  int use_qr = 0;
  get_blocked_qr(params, use_qr);
  elpa_set_integer(handle, "qr", use_qr, &error);
  
  /* Setup */
  assert_elpa_ok(elpa_setup(handle));

  /* Set tunables */
  constexpr int solver_enum = ELPA_SOLVER_2STAGE;
  elpa_set_integer(handle, "solver", solver_enum, &error);
  assert_elpa_ok(error);

  int kernel = get_kernel(params);
  elpa_set_integer(handle, "real_kernel", kernel, &error);
  assert_elpa_ok(error);

  // call eigenvalue routine
  int info = elpa::diag(handle, mat, eigvals);
  elpa_deallocate(handle, &error);

  params_out.set("info", info);
  return params_out;
}

} // namespace elpa
} // namespace rokko

#endif // ROKKO_ELPA_DIAGONALIZE_ELPA2_HPP
