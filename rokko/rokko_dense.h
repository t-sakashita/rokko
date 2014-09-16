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

#ifndef ROKKO_DENSE_H
#define ROKKO_DENSE_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  rokko_grid_col_major = 1, rokko_grid_row_major = 2,
  rokko_matrix_col_major = 3, rokko_matrix_row_major = 4
};

struct rokko_serial_dense_solver {
  void* ptr;
};

struct rokko_grid {
  void* ptr;
  int major;
};

struct rokko_parallel_dense_solver {
  void* ptr;
};

struct rokko_localized_vector {
  void* ptr;
};

struct rokko_localized_matrix {
  void* ptr;
  int major;
};

struct rokko_distributed_matrix {
  void* ptr;
  int major;
};


/* serial_dense_sover */
void rokko_serial_dense_solver_construct(struct rokko_serial_dense_solver* solver, char* solver_name, int argc, char** argv);
void rokko_serial_dense_solver_construct_f(struct rokko_serial_dense_solver* solver, char* solver_name);
void rokko_serial_dense_solver_destruct(struct rokko_serial_dense_solver* solver);
void rokko_serial_dense_solver_diagonalize_localized_matrix(struct rokko_serial_dense_solver* solver,
							    struct rokko_localized_matrix* mat, struct rokko_localized_vector* eigvals,
							    struct rokko_localized_matrix* eigvecs);


/* parallel_dense_sover */

void rokko_grid_construct(struct rokko_grid* grid, MPI_Comm comm, int grid_major);

void rokko_grid_construct_f(struct rokko_grid* grid, int comm, int grid_major);

void rokko_grid_destruct(struct rokko_grid* grid);

int rokko_grid_get_myrank(struct rokko_grid grid);

int rokko_grid_get_nprocs(struct rokko_grid grid);
 
void rokko_parallel_dense_solver_construct(struct rokko_parallel_dense_solver* solver, char* solver_name, int argc, char** argv);

void rokko_parallel_dense_solver_construct_f(struct rokko_parallel_dense_solver* solver, char* solver_name);

void rokko_parallel_dense_solver_destruct(struct rokko_parallel_dense_solver* solver);

void rokko_parallel_dense_solver_diagonalize_distributed_matrix(struct rokko_parallel_dense_solver* solver,
						 struct rokko_distributed_matrix* mat, struct rokko_localized_vector* eigvals,
						 struct rokko_distributed_matrix* eigvecs);

/* localized_matrix */

void rokko_localized_vector_construct(struct rokko_localized_vector* vec, int dim1);

void rokko_localized_vector_destruct(struct rokko_localized_vector* vec);

double rokko_localized_vector_get(struct rokko_localized_vector vec, int i);
  
double rokko_localized_vector_get_f(struct rokko_localized_vector vec, int i);

void rokko_localized_matrix_construct(struct rokko_localized_matrix* matrix, int dim1, int dim2,
				      int matrix_major);
  
void rokko_localized_matrix_destruct(struct rokko_localized_matrix* matrix);

void rokko_localized_matrix_print(struct rokko_localized_matrix matrix);

/* distributed_matrix */
void rokko_distributed_matrix_construct(struct rokko_distributed_matrix* matrix, int dim1, int dim2,
					struct rokko_grid grid, struct rokko_parallel_dense_solver solver, int matrix_major);
  
void rokko_distributed_matrix_destruct(struct rokko_distributed_matrix* matrix);

void rokko_distributed_matrix_generate_function(struct rokko_distributed_matrix* matrix,
						double (*func)(int i, int j));

void rokko_distributed_matrix_print(struct rokko_distributed_matrix matrix);

void rokko_distributed_matrix_set_local(struct rokko_distributed_matrix* matrix, int local_i, int local_j, double value);

double rokko_distributed_matrix_get_local(struct rokko_distributed_matrix matrix, int local_i, int local_j);

void rokko_distributed_matrix_set_global(struct rokko_distributed_matrix* matrix, int global_i, int global_j, double value);

double rokko_distributed_matrix_get_global(struct rokko_distributed_matrix matrix, int global_i, int global_j);

int rokko_distributed_matrix_get_m_local(struct rokko_distributed_matrix matrix);

int rokko_distributed_matrix_get_n_local(struct rokko_distributed_matrix matrix);

int rokko_distributed_matrix_get_m_global(struct rokko_distributed_matrix matrix);

int rokko_distributed_matrix_get_n_global(struct rokko_distributed_matrix matrix);

int rokko_distributed_matrix_get_nprocs(struct rokko_distributed_matrix matrix);

int rokko_distributed_matrix_get_myrank(struct rokko_distributed_matrix matrix);

int rokko_distributed_matrix_translate_l2g_row(struct rokko_distributed_matrix matrix, int local_i);

int rokko_distributed_matrix_translate_l2g_col(struct rokko_distributed_matrix matrix, int local_j);

int rokko_distributed_matrix_translate_g2l_row(struct rokko_distributed_matrix matrix, int global_i);

int rokko_distributed_matrix_translate_g2l_col(struct rokko_distributed_matrix matrix, int global_j);

void rokko_gather(struct rokko_distributed_matrix* matrix, double* array, int root);

void rokko_scatter(double* global_array, struct rokko_distributed_matrix* matrix, int root);

void rokko_all_gather(struct rokko_distributed_matrix* matrix, double* array);

#ifdef __cplusplus
}
#endif

#endif // ROKKO_DENSE_H
