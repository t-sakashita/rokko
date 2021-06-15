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

#pragma once

#include <rokko/eigen_matrix.h>
#include <rokko/distributed_matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

void rokko_gather(struct rokko_distributed_matrix matrix, double* array, int root);

void rokko_scatter(double* global_array, struct rokko_distributed_matrix matrix, int root);

void rokko_gather_eigen_matrix(struct rokko_distributed_matrix matrix, struct rokko_eigen_matrix lmatrix, int root);
void rokko_scatter_eigen_matrix(struct rokko_eigen_matrix lmatrix, struct rokko_distributed_matrix matrix, int root);

void rokko_all_gather(struct rokko_distributed_matrix matrix, double* array);

#ifdef __cplusplus
}
#endif
