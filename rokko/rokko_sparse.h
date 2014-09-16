/*****************************************************************************
 *
 * Rokko: Integrated Interface for libraries of eigenvalue decomposition
 *
 * Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
 *                            Synge Todo <wistaria@comp-phys.org>,
 *                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
 *
 * Distributed under the Boost Software License, Version 1.0. (See accompanying
 * file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *****************************************************************************/

#ifndef ROKKO_SPARSE_H
#define ROKKO_SPARSE_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct rokko_parallel_sparse_solver {
  void* ptr;
};

struct rokko_distributed_crs_matrix {
  void* ptr;
};

void rokko_parallel_sparse_solver_construct(struct rokko_parallel_sparse_solver* solver, char* solver_name, int argc, char** argv);

void rokko_parallel_sparse_solver_construct_f(struct rokko_parallel_sparse_solver* solver, char* solver_name);

void rokko_parallel_sparse_solver_destruct(struct rokko_parallel_sparse_solver* solver);

void rokko_parallel_sparse_solver_diagonalize_distributed_crs_matrix(struct rokko_parallel_sparse_solver* solver,
								 struct rokko_distributed_crs_matrix* mat,
								 int num_evals, int block_size, int max_iters, double tol);

double rokko_parallel_sparse_solver_eigenvalue(struct rokko_parallel_sparse_solver* solver, int i);

void rokko_parallel_sparse_solver_eigenvector(struct rokko_parallel_sparse_solver* solver, int i, double* vec);

int rokko_parallel_sparse_solver_num_conv(struct rokko_parallel_sparse_solver* solver);

/* rokko_distributed_crs_matrix */
void rokko_distributed_crs_matrix_construct(struct rokko_distributed_crs_matrix* matrix, int dim1, int dim2,
					    struct rokko_parallel_sparse_solver solver);
  
void rokko_distributed_crs_matrix_destruct(struct rokko_distributed_crs_matrix* matrix);

void rokko_distributed_crs_matrix_insert(struct rokko_distributed_crs_matrix* matrix, int row, int col_size, int* cols, double* values);

void rokko_distributed_crs_matrix_complete(struct rokko_distributed_crs_matrix* matrix);

int rokko_distributed_crs_matrix_num_local_rows(struct rokko_distributed_crs_matrix* matrix);

int rokko_distributed_crs_matrix_start_row(struct rokko_distributed_crs_matrix* matrix);

int rokko_distributed_crs_matrix_end_row(struct rokko_distributed_crs_matrix* matrix);

void rokko_distributed_crs_matrix_print(struct rokko_distributed_crs_matrix* matrix);


#ifdef __cplusplus
}
#endif

#endif // ROKKO_SPARSE_H
