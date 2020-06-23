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

#ifndef ROKKO_ELPA_DIAGONALIZE_ELPA1_HPP
#define ROKKO_ELPA_DIAGONALIZE_ELPA1_HPP

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
parameters diagonalize_elpa1(distributed_matrix<T, MATRIX_MAJOR>& mat,
			     VEC& eigvals, distributed_matrix<T, MATRIX_MAJOR>& eigvecs,
			     parameters const& params) {
  parameters params_out;
  if(mat.is_row_major())
    throw std::invalid_argument("elpa::diagonalize_elpa1() : elpa doesn't support matrix_row_major.  Use elpa with matrix_col_major.");

  int error;
  elpa_t handle = elpa_allocate(&error);

  set_parameters(mat, params, handle);

  /* Setup */
  assert_elpa_ok(elpa_setup(handle));

  /* Set tunables */
  constexpr int solver_enum = ELPA_SOLVER_1STAGE;
  elpa_set_integer(handle, "solver", solver_enum, &error);
  assert_elpa_ok(error);
  
  int info = elpa::diag(handle, mat, eigvals, eigvecs);
  elpa_deallocate(handle, &error);

  params_out.set("info", info);
  return params_out;
}

template<typename T, typename MATRIX_MAJOR, typename VEC>
parameters diagonalize_elpa1(distributed_matrix<T, MATRIX_MAJOR>& mat,
			     VEC& eigvals,
			     parameters const& params) {
  parameters params_out;
  if(mat.is_row_major())
    throw std::invalid_argument("elpa::diagonalize_elpa1() : elpa doesn't support matrix_row_major.  Use elpa with matrix_col_major.");

  int error;
  elpa_t handle = elpa_allocate(&error);

  /* Set parameters */
  set_parameters(mat, params, handle);
  
  /* Setup */
  assert_elpa_ok(elpa_setup(handle));

  /* Set tunables */
  constexpr int solver_enum = ELPA_SOLVER_1STAGE;
  elpa_set_integer(handle, "solver", solver_enum, &error);
  assert_elpa_ok(error);

  int info = elpa::diag(handle, mat, eigvals);
  elpa_deallocate(handle, &error);

  params_out.set("info", info);
  return params_out;
}

} // namespace elpa
} // namespace rokko

#endif // ROKKO_ELPA_DIAGONALIZE_ELPA1_HPP
