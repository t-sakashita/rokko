/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MAPPING_BC_H
#define ROKKO_MAPPING_BC_H

#include <rokko/grid.h>
#include <rokko/matrix_major.h>

#ifdef __cplusplus
extern "C" {
#endif

struct rokko_mapping_bc {
  void* ptr;
  int major;
};

void rokko_mapping_bc_construct_block_size(struct rokko_mapping_bc* map, int global_dim, int block_size, struct rokko_grid grid);
void rokko_mapping_bc_destruct(struct rokko_mapping_bc* map);
int rokko_mapping_bc_get_mb(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_nb(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_m_local(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_n_local(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_m_global(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_n_global(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_lld(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_m_size(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_n_size(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_length_array(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_nprocs(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_myrank(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_nprow(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_npcol(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_myrow(struct rokko_mapping_bc map);
int rokko_mapping_bc_get_mycol(struct rokko_mapping_bc map);
int rokko_mapping_bc_translate_l2g_row(struct rokko_mapping_bc map, int local_i);
int rokko_mapping_bc_translate_l2g_col(struct rokko_mapping_bc map, int local_j);
int rokko_mapping_bc_translate_g2l_row(struct rokko_mapping_bc map, int global_i);
int rokko_mapping_bc_translate_g2l_col(struct rokko_mapping_bc map, int global_j);
/* offset by one for Fortran */
int rokko_mapping_bc_translate_l2g_row1(struct rokko_mapping_bc map, int local_i);
int rokko_mapping_bc_translate_l2g_col1(struct rokko_mapping_bc map, int local_j);
int rokko_mapping_bc_translate_g2l_row1(struct rokko_mapping_bc map, int global_i);
int rokko_mapping_bc_translate_g2l_col1(struct rokko_mapping_bc map, int global_j);
bool rokko_mapping_bc_is_row_major(struct rokko_mapping_bc map);
bool rokko_mapping_bc_is_col_major(struct rokko_mapping_bc map);

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_MAPPING_BC_H */
