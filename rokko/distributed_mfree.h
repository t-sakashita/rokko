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

#ifndef ROKKO_DISTRIBUTED_MFREE_H
#define ROKKO_DISTRIBUTED_MFREE_H

#ifdef __cplusplus
extern "C" {
#endif

struct rokko_distributed_mfree {
  void* ptr;
};

void rokko_distributed_mfree_construct(struct rokko_distributed_mfree* matrix,
				       void (*multiply)(const double *const, double *const, void*),
				       void* vars,
				       int dim, MPI_Comm comm);
void rokko_distributed_mfree_construct_no_context(struct rokko_distributed_mfree* matrix,
				       void (*multiply)(const double *const, double *const),
				       int dim, MPI_Comm comm);
void rokko_distributed_mfree_destruct(struct rokko_distributed_mfree* matrix);
void rokko_distributed_mfree_f_construct(struct rokko_distributed_mfree* matrix,
					 void (*multiply)(const int*, const double *const, double *const),
					 int dim, int comm_f);
void rokko_distributed_mfree_f_destruct(struct rokko_distributed_mfree* matrix);
int rokko_distributed_mfree_dim(struct rokko_distributed_mfree matrix);
int rokko_distributed_mfree_num_local_rows(struct rokko_distributed_mfree matrix);
int rokko_distributed_mfree_start_row(struct rokko_distributed_mfree matrix);
int rokko_distributed_mfree_end_row(struct rokko_distributed_mfree matrix);
int rokko_distributed_mfree_offset(struct rokko_distributed_mfree matrix);

int rokko_distributed_mfree_f_dim(struct rokko_distributed_mfree* matrix);
int rokko_distributed_mfree_f_num_local_rows(struct rokko_distributed_mfree* matrix);
int rokko_distributed_mfree_f_start_row(struct rokko_distributed_mfree* matrix);
int rokko_distributed_mfree_f_end_row(struct rokko_distributed_mfree* matrix);
int rokko_distributed_mfree_f_offset(struct rokko_distributed_mfree* matrix);

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_DISTRIBUTED_MFREE_H */
