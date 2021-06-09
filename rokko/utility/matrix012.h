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

#pragma once

#include <rokko/dense.h>

#ifdef __cplusplus
extern "C"{
#endif

void rokko_matrix012_generate_eigen_matrix(struct rokko_eigen_matrix matrix);

#if defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
void rokko_matrix012_generate_distributed_matrix(struct rokko_distributed_matrix matrix);
#endif /* defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER) */

#ifdef __cplusplus
}
#endif
