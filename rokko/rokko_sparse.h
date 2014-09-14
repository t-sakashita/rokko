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

void rokko_paralle_sparse_solver_construct_f(struct rokko_parallel_sparse_solver* solver, char* solver_name);

void rokko_parallel_sparse_solver_destruct(struct rokko_parallel_sparse_solver* solver);

/*
void rokko_parallel_sparse_solver_diagonalize_distributed_crs_matrix(struct rokko_parallel_sparse_solver* solver,
								     struct rokko_distributed_crs_matrix* mat, struct rokko_localized_vector* eigvals,
								     struct rokko_distributed_crs_matrix* eigvecs);
*/

void rokko_distributed_crs_matrix_construct(struct rokko_distributed_crs_matrix* matrix, int dim1, int dim2,
					    struct rokko_parallel_sparse_solver solver);
  
void rokko_distributed_crs_matrix_destruct(struct rokko_distributed_crs_matrix* matrix);

void rokko_distributed_crs_matrix_insert(struct rokko_distributed_crs_matrix* matrix, int row, int col_size, int* cols, double* values);

int rokko_distributed_crs_matrix_start_row(struct rokko_distributed_crs_matrix* matrix);

int rokko_distributed_crs_matrix_end_row(struct rokko_distributed_crs_matrix* matrix);

void rokko_distributed_crs_matrix_complete(struct rokko_distributed_crs_matrix* matrix);

/*
void rokko_distributed_crs_matrix_generate_function(struct rokko_distributed_crs_matrix* matrix,
						    double (*func)(int i, int j));

void rokko_distributed_crs_matrix_print(struct rokko_distributed_crs_matrix matrix);

void rokko_distributed_crs_matrix_set_local(struct rokko_distributed_crs_matrix* matrix, int local_i, int local_j, double value);

double rokko_distributed_crs_matrix_get_local(struct rokko_distributed_crs_matrix matrix, int local_i, int local_j);

void rokko_distributed_crs_matrix_set_global(struct rokko_distributed_crs_matrix* matrix, int global_i, int global_j, double value);

double rokko_distributed_crs_matrix_get_global(struct rokko_distributed_crs_matrix matrix, int global_i, int global_j);

int rokko_distributed_crs_matrix_get_m_local(struct rokko_distributed_crs_matrix matrix);

int rokko_distributed_crs_matrix_get_n_local(struct rokko_distributed_crs_matrix matrix);

int rokko_distributed_crs_matrix_get_m_global(struct rokko_distributed_crs_matrix matrix);

int rokko_distributed_crs_matrix_get_n_global(struct rokko_distributed_crs_matrix matrix);

int rokko_distributed_crs_matrix_get_nprocs(struct rokko_distributed_crs_matrix matrix);

int rokko_distributed_crs_matrix_get_myrank(struct rokko_distributed_crs_matrix matrix);

int rokko_distributed_crs_matrix_translate_g2l_row(struct rokko_distributed_crs_matrix matrix, int global_i);

int rokko_distributed_crs_matrix_translate_g2l_col(struct rokko_distributed_crs_matrix matrix, int global_j);
*/

#ifdef __cplusplus
}
#endif

#endif // ROKKO_SPARSE_H
