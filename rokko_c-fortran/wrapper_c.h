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
namespace rokko {
extern "C"{
#endif

  void* initialize_distributed_matrix_row_major(int, int, void*, void*);
  void* initialize_distributed_matrix_col_major(int, int, void*, void*);
  void delete_distributed_matrix_col_major(void* );
  void delete_distributed_matrix_row_major(void* );

  void* initialize_localized_vector(int);
  void delete_localized_vector(void*);
  double localized_vector_get_element(void*, int);

  void* initialize_grid_col_major(MPI_Comm);
  void* initialize_grid_row_major(MPI_Comm);
  int grid_get_myrank(void*);
  int grid_get_nprocs(void*);
  void delete_grid(void*);

  void* initialize_solver(char* ,int, char*[]);
  void delete_solver(void* );
  void solver_diagonalize_matrix_col_major(void*, void*, void*, void*, void*);
  void solver_diagonalize_matrix_row_major(void*, void*, void*, void*, void*);

  void generate_distributed_matrix_function_col_major(void* mat, double (*func)(int i, int j));
  void generate_distributed_matrix_function_row_major(void* mat, double (*func)(int i, int j));

  void generate_distributed_matrix_function_col_major_fortran(void* mat, double (*func)(int* i, int* j));
  void generate_distributed_matrix_function_row_major_fortran(void* mat, double (*func)(int* i, int* j));

  void set_distributed_matrix_local_col_major(void* mat, int i, int j, double val);
  void set_distributed_matrix_local_row_major(void* mat, int i, int j, double val);

  double get_distributed_matrix_local_col_major(void* mat, int i, int j);
  double get_distributed_matrix_local_row_major(void* mat, int i, int j);

  void print_distributed_matrix_col_major(void* mat);
  void print_distributed_matrix_row_major(void* mat);

#ifdef __cplusplus
}
}
#endif
