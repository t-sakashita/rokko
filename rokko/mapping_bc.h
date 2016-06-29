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

#ifndef ROKKO_MAPPING_BC_H
#define ROKKO_MAPPING_BC_H

#include <rokko/parallel_dense_ev.h>

#ifdef __cplusplus
extern "C" {
#endif

/*void rokko_mapping_bc_construct(struct rokko_mapping_bc* map, int global_dim, struct rokko_grid grid, struct rokko_parallel_dense_ev solver);*/
void rokko_mapping_bc_construct_block_size(struct rokko_mapping_bc* map, int global_dim, int block_size, struct rokko_grid grid);
void rokko_mapping_bc_destruct(struct rokko_mapping_bc* map);

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_MAPPING_BC_H */
