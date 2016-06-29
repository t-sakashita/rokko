/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SERIAL_DENSE_EV_H
#define ROKKO_SERIAL_DENSE_EV_H

#ifdef __cplusplus
extern "C" {
#endif

struct rokko_serial_dense_ev {
  void* ptr;
};

void rokko_serial_dense_ev_construct(struct rokko_serial_dense_ev* solver,
  const char* solver_name, int argc, char** argv);
void rokko_serial_dense_ev_destruct(struct rokko_serial_dense_ev* solver);
struct rokko_parameters rokko_serial_dense_ev_diagonalize_localized_matrix(
  struct rokko_serial_dense_ev solver, struct rokko_localized_matrix mat,
  struct rokko_localized_vector eigval, struct rokko_localized_matrix eigvecs);
struct rokko_parameters rokko_serial_dense_ev_diagonalize_eigvals(
  struct rokko_serial_dense_ev solver, struct rokko_localized_matrix mat,
  struct rokko_localized_vector eigval);
struct rokko_parameters rokko_serial_dense_ev_diagonalize_eigvals_params(
  struct rokko_serial_dense_ev solver, struct rokko_localized_matrix mat,
  struct rokko_localized_vector eigvals, struct rokko_parameters params);
struct rokko_parameters rokko_serial_dense_ev_diagonalize_params(
  struct rokko_serial_dense_ev solver, struct rokko_localized_matrix mat,
  struct rokko_localized_vector eigvals, struct rokko_localized_matrix eigvecs,
  struct rokko_parameters params);

/* For Fortran binding */
void rokko_serial_dense_ev_construct_f(struct rokko_serial_dense_ev* solver,
  const char* solver_name);
void rokko_serial_dense_ev_diagonalize_f(struct rokko_serial_dense_ev* solver,
					 struct rokko_localized_matrix* mat, struct rokko_localized_vector* eigvals,
					 struct rokko_localized_matrix* eigvecs,
					 struct rokko_parameters* params, struct rokko_parameters* params_out);

void rokko_serial_dense_ev_diagonalize_no_params_out_f(struct rokko_serial_dense_ev* solver,
						       struct rokko_localized_matrix* mat, struct rokko_localized_vector* eigvals,
						       struct rokko_localized_matrix* eigvecs,
						       struct rokko_parameters* params);

void rokko_serial_dense_ev_diagonalize_no_params_inout_f(struct rokko_serial_dense_ev* solver,
						       struct rokko_localized_matrix* mat, struct rokko_localized_vector* eigvals,
						       struct rokko_localized_matrix* eigvecs);

void rokko_serial_dense_ev_diagonalize_eigvals_f(struct rokko_serial_dense_ev* solver,
						 struct rokko_localized_matrix* mat, struct rokko_localized_vector* eigvals,
						 struct rokko_parameters* params, struct rokko_parameters* params_out);

void rokko_serial_dense_ev_diagonalize_eigvals_no_params_out_f(struct rokko_serial_dense_ev* solver,
							       struct rokko_localized_matrix* mat, struct rokko_localized_vector* eigvals,
							       struct rokko_parameters* params);

void rokko_serial_dense_ev_diagonalize_eigvals_no_params_inout_f(struct rokko_serial_dense_ev* solver,
								 struct rokko_localized_matrix* mat, struct rokko_localized_vector* eigvals);


int rokko_serial_dense_ev_num_solvers();
char** rokko_serial_dense_ev_solvers();
char* rokko_serial_dense_ev_default_solver();

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_SERIAL_DENSE_EV_H */
