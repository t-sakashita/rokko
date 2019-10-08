/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include <rokko/pyrokko_distributed_mfree.hpp>
#include <rokko/utility/heisenberg_hamiltonian_mfree.hpp>

namespace py = pybind11;

PYBIND11_MODULE(heisenberg, m) {
  py::module::import("pyrokko");

  py::class_<rokko::heisenberg_mfree, rokko::distributed_mfree>(m, "mfree")
    .def(py::init<int, const std::vector<std::pair<int, int>>&>());
}