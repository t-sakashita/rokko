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

#include <rokko/parallel_dense_ev.hpp>
#include <rokko/eigenexa/solver.hpp>

ROKKO_REGISTER_PARALLEL_DENSE_SOLVER(rokko::eigenexa::solver, "eigenexa", 40)
