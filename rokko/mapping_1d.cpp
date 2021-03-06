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

#include "mapping_1d.hpp"

template<>
rokko::detail::ps_mapping_1d_factory *rokko::detail::ps_mapping_1d_factory::instance_ = nullptr;

template<>
rokko::detail::ps_mapping_1d_factory_num *rokko::detail::ps_mapping_1d_factory_num::instance_ = nullptr;
