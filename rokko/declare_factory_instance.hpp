/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2023 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <rokko/serial_dense_ev.hpp>
#include <rokko/parallel_dense_ev.hpp>
#include <rokko/parallel_sparse_ev.hpp>



#define ROKKO_DECLARE_FACTORY_INSTANCE(T) \
template<> \
std::unique_ptr<T>  T::instance_;


ROKKO_DECLARE_FACTORY_INSTANCE(rokko::detail::sd_solver_factory)

#ifdef ROKKO_HAVE_PARALLEL_DENSE_SOLVER
ROKKO_DECLARE_FACTORY_INSTANCE(rokko::detail::pd_solver_factory)
#endif

#ifdef ROKKO_HAVE_PARALLEL_SPARSE_SOLVER
ROKKO_DECLARE_FACTORY_INSTANCE(rokko::detail::ps_solver_factory)
ROKKO_DECLARE_FACTORY_INSTANCE(rokko::detail::ps_mapping_1d_factory)
ROKKO_DECLARE_FACTORY_INSTANCE(rokko::detail::ps_mapping_1d_factory_num)
ROKKO_DECLARE_FACTORY_INSTANCE(rokko::detail::ps_crs_factory)
#endif
