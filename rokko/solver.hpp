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

#ifndef ROKKO_SOLVER_HPP
#define ROKKO_SOLVER_HPP

#include <rokko/config.h>
#include <rokko/serial_dense_ev.hpp>
#ifdef ROKKO_HAVE_PARALLEL_DENSE_SOLVER
# include <rokko/parallel_dense_solver.hpp>
#endif
#ifdef ROKKO_HAVE_PARALLEL_SPARSE_SOLVER
# include <rokko/parallel_sparse_solver.hpp>
#endif

#endif // ROKKO_SOLVER_HPP
