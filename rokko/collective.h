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

#ifndef ROKKO_COLLECTIVE_H
#define ROKKO_COLLECTIVE_H

#ifdef __cplusplus
extern "C" {
#endif

void rokko_gather(struct rokko_distributed_matrix matrix, double* array, int root);

void rokko_scatter(double* global_array, struct rokko_distributed_matrix matrix, int root);

void rokko_gather_localized_matrix(struct rokko_distributed_matrix matrix, struct rokko_localized_matrix lmatrix, int root);
void rokko_scatter_localized_matrix(struct rokko_localized_matrix lmatrix, struct rokko_distributed_matrix matrix, int root);

void rokko_all_gather(struct rokko_distributed_matrix matrix, double* array);

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_COLLECTIVE_H */
