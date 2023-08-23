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

#pragma once

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
parameters diagonalize(distributed_matrix<T, rokko::matrix_col_major>& mat,
			     VEC& eigvals, distributed_matrix<T, rokko::matrix_col_major>& eigvecs,
			     parameters const& params) {
  parameters params_out;

  int error;
  const elpa_t handle = elpa_allocate(&error);

  set_parameters(mat, params, handle);
  assert_elpa_ok(elpa_setup(handle));

  set_solver(params, handle);
  const auto info = elpa::diag(handle, mat, eigvals, eigvecs);
  assert_elpa_ok( deallocate(handle) );

  params_out.set("info", info);
  return params_out;
}

template<typename T, typename VEC>
parameters diagonalize(distributed_matrix<T, rokko::matrix_row_major>& mat,
			     VEC& eigvals, distributed_matrix<T, rokko::matrix_row_major>& eigvecs,
			     parameters const& params) {
  throw std::invalid_argument("elpa::diagonalize_elpa2() : elpa doesn't support matrix_row_major.  Use elpa with matrix_col_major.");
}

template<typename T, typename VEC>
parameters diagonalize(distributed_matrix<T, rokko::matrix_col_major>& mat,
			     VEC& eigvals,
			     parameters const& params) {
  parameters params_out;

  int error;
  const elpa_t handle = elpa_allocate(&error);

  set_parameters(mat, params, handle);
  assert_elpa_ok(elpa_setup(handle));

  set_solver(params, handle);
  const auto info = elpa::diag(handle, mat, eigvals);
  assert_elpa_ok( deallocate(handle) );

  params_out.set("info", info);
  return params_out;
}

template<typename T, typename VEC>
parameters diagonalize(distributed_matrix<T, rokko::matrix_row_major>& mat,
			     VEC& eigvals,
			     parameters const& params) {
  throw std::invalid_argument("elpa::diagonalize_elpa2() : elpa doesn't support matrix_row_major.  Use elpa with matrix_col_major.");
}

} // namespace elpa
} // namespace rokko
