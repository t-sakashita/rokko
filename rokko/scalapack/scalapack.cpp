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

#include <rokko/parallel_dense_solver.hpp>
#include <rokko/scalapack/core.hpp>

ROKKO_REGISTER_PARALLEL_DENSE_SOLVER(rokko::scalapack::solver, "scalapack", 20)
