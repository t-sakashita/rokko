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

#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <rokko/pyrokko_communicator.hpp>
#include <rokko/grid.hpp>
#include <rokko/utility/tuple_to_array.hpp>

#include <memory>

namespace rokko {

namespace py = pybind11;

class wrap_grid : public grid {
public:
  using grid::grid;

  wrap_grid() : grid(to_MPI_Comm()) {}

  template <typename GRID_MAJOR>
  wrap_grid(GRID_MAJOR const& grid_major = grid_row_major) : grid(to_MPI_Comm(), grid_major) {}

  template <typename GRID_MAJOR>
  wrap_grid(pybind11::handle const& comm_handle, std::tuple<int,int> const& size_in, GRID_MAJOR const& grid_major = grid_row_major)
    : grid(to_MPI_Comm(comm_handle), to_array(size_in), grid_major) {}

  template <typename GRID_MAJOR>
  wrap_grid(pybind11::handle const& comm_handle, GRID_MAJOR const& grid_major = grid_row_major)
    : grid(to_MPI_Comm(comm_handle), grid_major) {}

  wrap_grid(rokko::grid const& g) : grid(g) {}

  std::string get_major_string() const {
    return grid::is_row_major() ? "row" : "col";
  }

  pybind11::handle get_comm() const {
    return pybind11::handle(PyMPIComm_New(grid::get_comm()));
  }

  std::tuple<int,int> get_shape() const {
    return grid::get_size();
  }

  std::tuple<int,int> get_mine() const {
    return grid::get_my_coordinate();
  }
};

} // end namespace rokko
