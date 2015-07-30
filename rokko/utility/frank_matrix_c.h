/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_FRANK_MATRIX_C_H
#define ROKKO_FRANK_MATRIX_C_H

#include <rokko/rokko_dense.h>

#ifdef __cplusplus
extern "C"{
#endif

void rokko_frank_matrix_generate_localized_matrix(struct rokko_localized_matrix* matrix);

#if defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
void rokko_frank_matrix_generate_distributed_matrix(struct rokko_distributed_matrix* matrix);
#endif


#ifdef __cplusplus
}
#endif

#endif
