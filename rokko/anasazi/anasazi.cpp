/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/parallel_sparse_ev.hpp>
#include <rokko/anasazi/solver.hpp>
#include <rokko/anasazi/mapping_1d.hpp>
#include <rokko/anasazi/distributed_crs_matrix.hpp>

ROKKO_REGISTER_PARALLEL_SPARSE_SOLVER(rokko::anasazi::solver, "anasazi", 40)

ROKKO_REGISTER_PARALLEL_SPARSE_MAPPING_1D(rokko::anasazi::mapping_1d, "anasazi", 40)

ROKKO_REGISTER_PARALLEL_SPARSE_CRS(rokko::anasazi::distributed_crs_matrix, "anasazi", 40)
