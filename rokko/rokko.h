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

#include <rokko/config.h>

#ifdef __cplusplus

#include <rokko/solver.hpp>
#if defined(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)
# include <rokko/grid_1d.hpp>
# include <rokko/distributed_vector.hpp>
#endif /* defined(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER) */

#else /* __cplus */

#include <rokko/solver_name.h>
#include <rokko/dense.h>
#if defined(ROKKO_HAVE_PARALLEL_SPARSE_SOLVER)
# include <rokko/grid.h>
# include <rokko/sparse.h>
#endif

#endif /* __cplus */
