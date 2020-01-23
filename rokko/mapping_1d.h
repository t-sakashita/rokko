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

#ifndef ROKKO_MAPPING_1D_H
#define ROKKO_MAPPING_1D_H

#ifdef __cplusplus
extern "C" {
#endif

struct rokko_mapping_1d {
  void* ptr;
};

void rokko_mapping_1d_construct(struct rokko_mapping_1d* map, int dim, MPI_Comm comm);
void rokko_mapping_1d_construct_f(struct rokko_mapping_1d* map, int dim, int comm_f);
void rokko_mapping_1d_destruct(struct rokko_mapping_1d* map);
int rokko_mapping_1d_get_dim(struct rokko_mapping_1d map);
int rokko_mapping_1d_num_local_rows(struct rokko_mapping_1d map);
int rokko_mapping_1d_start_row(struct rokko_mapping_1d map);
int rokko_mapping_1d_end_row(struct rokko_mapping_1d map);
MPI_Comm rokko_mapping_1d_get_comm(struct rokko_mapping_1d map);

MPI_Fint rokko_mapping_1d_get_comm_f(struct rokko_mapping_1d map);

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_MAPPING_1D_H */
