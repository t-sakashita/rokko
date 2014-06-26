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

#ifndef ROKKO_LAPACK_CPP
#define ROKKO_LAPACK_CPP

#include <rokko/solver_factory.hpp>
#include <rokko/lapack/core.hpp>

ROKKO_REGISTER_SERIAL_DENSE_SOLVER(rokko::lapack::solver<rokko::lapack::dsyev>, "lapack", 20)

#endif
