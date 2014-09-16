/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SOLVER_HPP
#define ROKKO_SOLVER_HPP

#include <rokko/config.hpp>

#include <rokko/serial_dense_solver.hpp>
#include <rokko/parallel_dense_solver.hpp>

#if defined(ROKKO_HAVE_ANASAZI) || defined(ROKKO_HAVE_SLEPC)
#include <rokko/parallel_sparse_solver.hpp>
#endif // defined(ROKKO_HAVE_ANASAZI) || defined(ROKKO_HAVE_SLEPC)

#endif // ROKKO_SOLVER_HPP
