/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Synge Todo <wistaria@comp-phys.org>,
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_EIGEN3_CPP
#define ROKKO_EIGEN3_CPP

#include <rokko/serial_dense_solver.hpp>
#include <rokko/eigen3/core.hpp>

ROKKO_REGISTER_SERIAL_DENSE_SOLVER(rokko::eigen3::solver, "eigen3", 10)

#endif
