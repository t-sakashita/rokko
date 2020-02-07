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

#ifndef ROKKO_SPARSE_H
#define ROKKO_SPARSE_H

#include <rokko/config.h>

#if defined(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)
# include <mpi.h>
#include <rokko/parallel_sparse_ev.h>
#include <rokko/mapping_1d.h>
#include <rokko/distributed_crs_matrix.h>
#include <rokko/distributed_mfree.h>

#endif /* defined(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER) */

#endif /* ROKKO_SPARSE_H */
