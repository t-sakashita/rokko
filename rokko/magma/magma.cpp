/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2014 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MAGMA_CPP
#define ROKKO_MAGMA_CPP

#include <rokko/serial_solver_factory.hpp>

#include <rokko/magma/core.hpp>

ROKKO_REGISTER_SERIAL_SOLVER(rokko::magma::solver<rokko::magma::dsyev>, "magma")

#endif // ROKKO_MAGMA_CPP
