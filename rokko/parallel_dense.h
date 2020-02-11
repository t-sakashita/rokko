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

#ifndef ROKKO_PARALLEL_DENSE_H
#define ROKKO_PARALLEL_DENSE_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <rokko/grid.h>
#include <rokko/mapping_bc.h>
#include <rokko/distributed_matrix.h>
#include <rokko/parallel_dense_ev.h>

#include <rokko/collective.h>

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_PARALLEL_DENSE_H */
