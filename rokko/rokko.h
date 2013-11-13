/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>,
*                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  rokko_grid_col_major = 1, rokko_grid_row_major = 2,
  rokko_matrix_col_major = 3, rokko_matrix_row_major = 4
};
  
struct rokko_grid {
  void* ptr;
  int major;
};

struct rokko_solver {
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

void rokko_grid_construct(struct rokko_grid* grid, MPI_Comm comm, int grid_major);

void rokko_grid_construct_f(struct rokko_grid* grid, int comm, int grid_major);

void rokko_grid_destruct(struct rokko_grid* grid);

int rokko_grid_get_myrank(struct rokko_grid grid);

int rokko_grid_get_nprocs(struct rokko_grid grid);
 
void rokko_solver_construct(struct rokko_solver* solver, char* solver_name, int argc, char** argv);

void rokko_solver_construct_f(struct rokko_solver* solver, char* solver_name);

void rokko_solver_destruct(struct rokko_solver* solver);

void rokko_solver_diagonalize_distributed_matrix(struct rokko_solver* solver,
  struct rokko_distributed_matrix* mat, struct rokko_localized_vector* eigvals,
  struct rokko_distributed_matrix* eigvecs);

void rokko_localized_vector_construct(struct rokko_localized_vector* vec, int dim1);

void rokko_localized_vector_destruct(struct rokko_localized_vector* vec);

double rokko_localized_vector_get(struct rokko_localized_vector vec, int i);
  
double rokko_localized_vector_get_f(struct rokko_localized_vector vec, int i);

void rokko_localized_matrix_construct(struct rokko_localized_matrix* matrix, int dim1, int dim2,
  int matrix_major);
  
void rokko_localized_matrix_destruct(struct rokko_localized_matrix* matrix);

void rokko_distributed_matrix_construct(struct rokko_distributed_matrix* matrix, int dim1, int dim2,
    struct rokko_grid grid, struct rokko_solver solver, int matrix_major);
  
void rokko_distributed_matrix_destruct(struct rokko_distributed_matrix* matrix);

void rokko_distributed_matrix_generate_function(struct rokko_distributed_matrix* matrix,
  double (*func)(int i, int j));

void rokko_distributed_matrix_print(struct rokko_distributed_matrix matrix);

void rokko_distributed_matrix_set_local(struct rokko_distributed_matrix* matrix, int local_i, int local_j, double value);

double rokko_distributed_matrix_get_local(struct rokko_distributed_matrix matrix, int local_i, int local_j);

void rokko_distributed_matrix_set_global(struct rokko_distributed_matrix* matrix, int global_i, int global_j, double value);

double rokko_distributed_matrix_get_global(struct rokko_distributed_matrix matrix, int global_i, int global_j);



#ifdef __cplusplus
}
#endif
