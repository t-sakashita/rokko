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

#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <rokko/pyrokko_communicator.hpp>
#include <rokko/mapping_1d.hpp>

namespace rokko {

class wrap_mapping_1d : public mapping_1d {
public:
  using mapping_1d::mapping_1d;
  wrap_mapping_1d(mapping_1d const& map) : mapping_1d(map) {}

  explicit wrap_mapping_1d(int dim, pybind11::handle const& comm_handle, std::string const& solver_name)
    : mapping_1d(dim, wrap_communicator{comm_handle}, solver_name) {}

  explicit wrap_mapping_1d(int dim, pybind11::handle const& comm_handle)
    : mapping_1d(dim, wrap_communicator{comm_handle}) {}

  ~wrap_mapping_1d() = default;

};

} // end namespace rokko
