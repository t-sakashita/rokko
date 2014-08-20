/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/parallel_sparse_solver.hpp>
#include <rokko/anasazi/core.hpp>

ROKKO_REGISTER_PARALLEL_SPARSE_SOLVER(rokko::anasazi::solver, "anasazi", 40)

