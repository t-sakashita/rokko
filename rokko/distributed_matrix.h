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

#ifndef ROKKO_DISTRIBUTED_MATRIX_H
#define ROKKO_DISTRIBUTED_MATRIX_H

#include <rokko/mapping_bc.h>

#ifdef __cplusplus
extern "C" {
#endif

struct rokko_distributed_matrix {
  void* ptr;
  int major;
};

void rokko_distributed_matrix_construct(struct rokko_distributed_matrix* matrix, struct rokko_mapping_bc map);
void rokko_distributed_matrix_construct_array(struct rokko_distributed_matrix* matrix, struct rokko_mapping_bc map, double *array);
void rokko_distributed_matrix_construct_array_sizes(struct rokko_distributed_matrix* matrix, struct rokko_mapping_bc map, int dim1, int dim2, double *array);
void rokko_distributed_matrix_destruct(struct rokko_distributed_matrix* matrix);
void rokko_distributed_matrix_generate_function(struct rokko_distributed_matrix matrix,
  double (*func)(int i, int j));
void rokko_distributed_matrix_generate_function_p(struct rokko_distributed_matrix matrix,
  double (*func)(const int* i, const int* j));
void rokko_distributed_matrix_generate_function1(struct rokko_distributed_matrix matrix,
  double (*func)(int i, int j));
void rokko_distributed_matrix_generate_function1_p(struct rokko_distributed_matrix matrix,
  double (*func)(const int* i, const int* j));
void rokko_distributed_matrix_print(struct rokko_distributed_matrix matrix);
void rokko_distributed_matrix_set_local(struct rokko_distributed_matrix matrix,
  int local_i, int local_j, double value);
double rokko_distributed_matrix_get_local(struct rokko_distributed_matrix matrix,
  int local_i, int local_j);
void rokko_distributed_matrix_set_global(struct rokko_distributed_matrix matrix,
  int global_i, int global_j, double value);
double rokko_distributed_matrix_get_global(struct rokko_distributed_matrix matrix,
  int global_i, int global_j);
/* offset by one for Fortran */
void rokko_distributed_matrix_set_local1(struct rokko_distributed_matrix matrix,
  int local_i, int local_j, double value);
double rokko_distributed_matrix_get_local1(struct rokko_distributed_matrix matrix,
  int local_i, int local_j);
void rokko_distributed_matrix_set_global1(struct rokko_distributed_matrix matrix,
  int global_i, int global_j, double value);
double rokko_distributed_matrix_get_global1(struct rokko_distributed_matrix matrix,
  int global_i, int global_j);
int rokko_distributed_matrix_get_mb(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_get_nb(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_get_m_local(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_get_n_local(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_get_m_global(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_get_n_global(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_get_lld(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_get_m_size(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_get_n_size(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_get_nprocs(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_get_myrank(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_get_nprow(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_get_npcol(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_get_myrow(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_get_mycol(struct rokko_distributed_matrix matrix);
int rokko_distributed_matrix_translate_l2g_row(struct rokko_distributed_matrix matrix,
  int local_i);
int rokko_distributed_matrix_translate_l2g_col(struct rokko_distributed_matrix matrix,
  int local_j);
int rokko_distributed_matrix_translate_g2l_row(struct rokko_distributed_matrix matrix,
  int global_i);
int rokko_distributed_matrix_translate_g2l_col(struct rokko_distributed_matrix matrix,
  int global_j);
/* offset by one for Fortran */
int rokko_distributed_matrix_translate_l2g_row1(struct rokko_distributed_matrix matrix,
  int local_i);
int rokko_distributed_matrix_translate_l2g_col1(struct rokko_distributed_matrix matrix,
  int local_j);
int rokko_distributed_matrix_translate_g2l_row1(struct rokko_distributed_matrix matrix,
  int global_i);
int rokko_distributed_matrix_translate_g2l_col1(struct rokko_distributed_matrix matrix,
  int global_j);
double* rokko_distributed_matrix_get_array_pointer(struct rokko_distributed_matrix matrix);

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_DISTRIBUTED_MATRIX_H */
