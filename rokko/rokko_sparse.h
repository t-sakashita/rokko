/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SPARSE_H
#define ROKKO_SPARSE_H

#include <rokko/config.h>
#include <rokko/parameters.h>

#if defined(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)
# include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#if defined(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)
struct rokko_distributed_crs_matrix {
  void* ptr;
};

struct rokko_distributed_mfree {
  void* ptr;
};
#endif

#ifdef ROKKO_HAVE_PARALLEL_SPARSE_SOLVER
struct rokko_parallel_sparse_ev {
  void* ptr;
};
#endif
  
#if defined(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)

/* rokko_distributed_crs_matrix */
void rokko_distributed_crs_matrix_construct(struct rokko_distributed_crs_matrix* matrix,
					    int dim1, int dim2, struct rokko_parallel_sparse_ev solver);
void rokko_distributed_crs_matrix_destruct(struct rokko_distributed_crs_matrix* matrix);
void rokko_distributed_crs_matrix_insert(struct rokko_distributed_crs_matrix* matrix, int row,
					 int col_size, int* cols, double* values);
void rokko_distributed_crs_matrix_complete(struct rokko_distributed_crs_matrix* matrix);
int rokko_distributed_crs_matrix_num_local_rows(struct rokko_distributed_crs_matrix* matrix);
int rokko_distributed_crs_matrix_start_row(struct rokko_distributed_crs_matrix* matrix);
int rokko_distributed_crs_matrix_end_row(struct rokko_distributed_crs_matrix* matrix);
void rokko_distributed_crs_matrix_print(struct rokko_distributed_crs_matrix* matrix);

/* rokko_distributed_mfree */
void rokko_distributed_mfree_construct(struct rokko_distributed_mfree* matrix,
				       void (*multiply)(const double*, double*, void*),
				       void* vars,
				       int dim, int num_local_rows);
void rokko_distributed_mfree_destruct(struct rokko_distributed_mfree* matrix);
void rokko_distributed_mfree_f_construct(struct rokko_distributed_mfree* matrix,
					 void (*multiply)(int, const double*, double*),
					 int dim, int num_local_rows);
void rokko_distributed_mfree_f_destruct(struct rokko_distributed_mfree* matrix);
int rokko_distributed_mfree_dim(struct rokko_distributed_mfree* matrix);
int rokko_distributed_mfree_num_local_rows(struct rokko_distributed_mfree* matrix);
int rokko_distributed_mfree_offset(struct rokko_distributed_mfree* matrix);

/* rokko_parallel_sparse_ev */
void rokko_parallel_sparse_ev_construct(struct rokko_parallel_sparse_ev* solver,
  const char* solver_name, int argc, char** argv);
void rokko_parallel_sparse_ev_construct_f(struct rokko_parallel_sparse_ev* solver,
  const char* solver_name);
void rokko_parallel_sparse_ev_destruct(struct rokko_parallel_sparse_ev* solver);
struct rokko_parameters rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix_noparams(
  struct rokko_parallel_sparse_ev* solver, struct rokko_distributed_crs_matrix* mat,
  int num_evals, int block_size, int max_iters, double tol);
struct rokko_parameters rokko_parallel_sparse_ev_diagonalize_distributed_mfree_noparams(
  struct rokko_parallel_sparse_ev* solver, struct rokko_distributed_mfree* mat,
  int num_evals, int block_size, int max_iters, double tol);
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

double rokko_parallel_sparse_ev_eigenvalue(struct rokko_parallel_sparse_ev* solver, int i);
void rokko_parallel_sparse_ev_eigenvector(struct rokko_parallel_sparse_ev* solver, int i,
  double* vec);
int rokko_parallel_sparse_ev_num_conv(struct rokko_parallel_sparse_ev* solver);

int rokko_parallel_sparse_ev_num_solvers();
char** rokko_parallel_sparse_ev_solvers();
char* rokko_parallel_sparse_ev_default_solver();

#endif

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_SPARSE_H */
