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

#ifndef ROKKO_DISTRIBUTED_CRS_MATRIX_H
#define ROKKO_DISTRIBUTED_CRS_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

struct rokko_distributed_crs_matrix {
  void* ptr;
};

struct rokko_parallel_sparse_ev;

void rokko_distributed_crs_matrix_construct_new(struct rokko_distributed_crs_matrix* matrix, struct rokko_mapping_1d map, int num_entries_per_row);
void rokko_distributed_crs_matrix_construct(struct rokko_distributed_crs_matrix* matrix,
					    int dim1, int dim2, struct rokko_parallel_sparse_ev solver);
void rokko_distributed_crs_matrix_destruct(struct rokko_distributed_crs_matrix* matrix);
void rokko_distributed_crs_matrix_insert(struct rokko_distributed_crs_matrix matrix, int row,
					 int col_size, int* cols, double* values);
void rokko_distributed_crs_matrix_complete(struct rokko_distributed_crs_matrix matrix);
int rokko_distributed_crs_matrix_num_local_rows(struct rokko_distributed_crs_matrix matrix);
int rokko_distributed_crs_matrix_start_row(struct rokko_distributed_crs_matrix matrix);
int rokko_distributed_crs_matrix_end_row(struct rokko_distributed_crs_matrix matrix);
void rokko_distributed_crs_matrix_print(struct rokko_distributed_crs_matrix matrix);

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_DISTRIBUTED_CRS_MATRIX_H */
