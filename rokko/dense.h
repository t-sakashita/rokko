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

#ifndef ROKKO_DENSE_H
#define ROKKO_DENSE_H

#include <rokko/config.h>

#include <rokko/serial_dense.h>

#if defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
#include <rokko/parallel_dense.h>
#endif /* defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER) */

#endif /* ROKKO_DENSE_H */
