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

#ifndef ROKKO_PARALLEL_SPARSE_EV_H
#define ROKKO_PARALLEL_SPARSE_EV_H

#include <rokko/mapping_1d.h>
#include <rokko/distributed_crs_matrix.h>
#include <rokko/distributed_mfree.h>
#include <rokko/parameters.h>

#ifdef __cplusplus
extern "C" {
#endif

struct rokko_parallel_sparse_ev {
  void* ptr;
};

void rokko_parallel_sparse_ev_construct(struct rokko_parallel_sparse_ev* solver,
  const char* solver_name, int argc, char** argv);
void rokko_parallel_sparse_ev_construct_f(struct rokko_parallel_sparse_ev* solver,
  const char* solver_name);
void rokko_parallel_sparse_ev_destruct(struct rokko_parallel_sparse_ev* solver);
struct rokko_parameters rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix(struct rokko_parallel_sparse_ev solver,
										    struct rokko_distributed_crs_matrix mat,
										    struct rokko_parameters params);
struct rokko_parameters rokko_parallel_sparse_ev_diagonalize_distributed_mfree(struct rokko_parallel_sparse_ev solver,
									       struct rokko_distributed_mfree mat,
									       struct rokko_parameters params);

/* Fortran binding */
void rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix_f(struct rokko_parallel_sparse_ev* solver,
								   struct rokko_distributed_crs_matrix* mat,
								   struct rokko_parameters* params, struct rokko_parameters* params_out);
void rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix_noreturn_f(struct rokko_parallel_sparse_ev* solver,
									    struct rokko_distributed_crs_matrix* mat,
									    struct rokko_parameters* params);
void rokko_parallel_sparse_ev_diagonalize_distributed_mfree_f(struct rokko_parallel_sparse_ev* solver,
							      struct rokko_distributed_mfree* mat,
							      struct rokko_parameters* params, struct rokko_parameters* params_out);
void rokko_parallel_sparse_ev_diagonalize_distributed_mfree_noreturn_f(struct rokko_parallel_sparse_ev* solver,
								       struct rokko_distributed_mfree* mat,
								       struct rokko_parameters* params);

double rokko_parallel_sparse_ev_eigenvalue(struct rokko_parallel_sparse_ev solver, int i);
void rokko_parallel_sparse_ev_eigenvector(struct rokko_parallel_sparse_ev solver, int i,
  double* vec);
int rokko_parallel_sparse_ev_num_conv(struct rokko_parallel_sparse_ev solver);

struct rokko_mapping_1d rokko_parallel_sparse_ev_default_mapping(struct rokko_parallel_sparse_ev solver, int dim, MPI_Comm comm);
void rokko_parallel_sparse_ev_default_mapping_f(struct rokko_parallel_sparse_ev solver, int dim, int comm_f, struct rokko_mapping_1d* map);

int rokko_parallel_sparse_ev_num_solvers();
char** rokko_parallel_sparse_ev_solvers();
char* rokko_parallel_sparse_ev_default_solver();

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_PARALLEL_SPARSE_EV_H */
