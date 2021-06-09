/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <rokko/grid.h>
#include <rokko/mapping_bc.h>
#include <rokko/distributed_matrix.h>
#include <rokko/eigen_vector.h>
#include <rokko/parameters.h>

#ifdef __cplusplus
extern "C" {
#endif

struct rokko_parallel_dense_ev {
  void* ptr;
};

void rokko_parallel_dense_ev_construct(struct rokko_parallel_dense_ev* solver,
  const char* solver_name, int argc, char** argv);
void rokko_parallel_dense_ev_construct_f(struct rokko_parallel_dense_ev* solver,
  const char* solver_name);
void rokko_parallel_dense_ev_destruct(struct rokko_parallel_dense_ev* solver);

struct rokko_parameters rokko_parallel_dense_ev_diagonalize(
  struct rokko_parallel_dense_ev solver, struct rokko_distributed_matrix mat,
  struct rokko_eigen_vector eigval, struct rokko_distributed_matrix eigvecs, struct rokko_parameters);
struct rokko_parameters rokko_parallel_dense_ev_diagonalize_distributed_matrix(
  struct rokko_parallel_dense_ev solver, struct rokko_distributed_matrix mat,
  struct rokko_eigen_vector eigval, struct rokko_distributed_matrix eigvecs);
struct rokko_parameters rokko_parallel_dense_ev_diagonalize_eigvals(
  struct rokko_parallel_dense_ev solver, struct rokko_distributed_matrix mat,
  struct rokko_eigen_vector eigval);
struct rokko_parameters rokko_parallel_dense_ev_diagonalize_eigvals_params(
  struct rokko_parallel_dense_ev solver, struct rokko_distributed_matrix mat,
  struct rokko_eigen_vector eigvals, struct rokko_parameters params);
struct rokko_parameters rokko_parallel_dense_ev_diagonalize_params(
  struct rokko_parallel_dense_ev solver, struct rokko_distributed_matrix mat,
  struct rokko_eigen_vector eigvals, struct rokko_distributed_matrix eigvecs,
  struct rokko_parameters params);

/* For Fortran binding */
void rokko_parallel_dense_ev_construct_f(struct rokko_parallel_dense_ev* solver,
  const char* solver_name);
void rokko_parallel_dense_ev_diagonalize_f(struct rokko_parallel_dense_ev* solver,
					 struct rokko_distributed_matrix* mat, struct rokko_eigen_vector* eigvals,
					 struct rokko_distributed_matrix* eigvecs,
					 struct rokko_parameters* params, struct rokko_parameters* params_out);

void rokko_parallel_dense_ev_diagonalize_no_params_out_f(struct rokko_parallel_dense_ev* solver,
						       struct rokko_distributed_matrix* mat, struct rokko_eigen_vector* eigvals,
						       struct rokko_distributed_matrix* eigvecs,
						       struct rokko_parameters* params);

void rokko_parallel_dense_ev_diagonalize_no_params_inout_f(struct rokko_parallel_dense_ev* solver,
						       struct rokko_distributed_matrix* mat, struct rokko_eigen_vector* eigvals,
						       struct rokko_distributed_matrix* eigvecs);

void rokko_parallel_dense_ev_diagonalize_eigvals_f(struct rokko_parallel_dense_ev* solver,
						 struct rokko_distributed_matrix* mat, struct rokko_eigen_vector* eigvals,
						 struct rokko_parameters* params, struct rokko_parameters* params_out);

void rokko_parallel_dense_ev_diagonalize_eigvals_no_params_out_f(struct rokko_parallel_dense_ev* solver,
							       struct rokko_distributed_matrix* mat, struct rokko_eigen_vector* eigvals,
							       struct rokko_parameters* params);

void rokko_parallel_dense_ev_diagonalize_eigvals_no_params_inout_f(struct rokko_parallel_dense_ev* solver,
								 struct rokko_distributed_matrix* mat, struct rokko_eigen_vector* eigvals);

struct rokko_mapping_bc rokko_parallel_dense_ev_default_mapping(struct rokko_parallel_dense_ev solver, int dim, struct rokko_grid g);
void rokko_parallel_dense_ev_default_mapping_f(struct rokko_parallel_dense_ev solver, int dim, struct rokko_grid g, struct rokko_mapping_bc* map);

int rokko_parallel_dense_ev_num_solvers();
char** rokko_parallel_dense_ev_solvers();
const char* rokko_parallel_dense_ev_default_solver();

#ifdef __cplusplus
}
#endif
