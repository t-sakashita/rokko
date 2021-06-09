/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include <rokko/utility/solver_name.hpp>

namespace rokko {

namespace py = pybind11;

std::tuple<std::string,std::string> wrap_split_solver_name(std::string const& str) {
  std::string library, routine;
  split_solver_name(str, library, routine);
  return std::tuple<std::string,std::string>(library, routine);
}

} // end namespace rokko
