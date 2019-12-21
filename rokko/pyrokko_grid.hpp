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

#ifndef PYROKKO_GRID_HPP
#define PYROKKO_GRID_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <rokko/pyrokko_communicator.hpp>
#include <rokko/grid.hpp>

namespace rokko {

namespace py = pybind11;

class wrap_grid {
public:
  wrap_grid() {
    const int rc = import_mpi4py();
    assert(rc==0);
    wrap_communicator wrap_comm;
    MPI_Comm comm = wrap_comm.get_comm();  // default: grid_row_major
    ptr = new grid(comm);
  }

  template <typename GRID_MAJOR>
  wrap_grid(GRID_MAJOR const& grid_major = grid_row_major) {
    const int rc = import_mpi4py();
    assert(rc==0);
    wrap_communicator wrap_comm;
    MPI_Comm comm = wrap_comm.get_comm();
    ptr = new grid(comm, grid_major);
  }
  
  template <typename GRID_MAJOR>
  wrap_grid(pybind11::handle const& comm_handle, GRID_MAJOR const& grid_major = grid_row_major, int lld = 0) {
    const int rc = import_mpi4py();
    assert(rc==0);
    wrap_communicator wrap_comm(*PyMPIComm_Get(comm_handle.ptr()));
    MPI_Comm comm = wrap_comm.get_comm();
    ptr = new grid(comm, grid_major, lld);
  }

  grid const& get_grid() const {
    return *ptr;
  }

  std::string get_major_string() const {
    if (get_grid().is_row_major())  return "row";
    else  return "col";
  }

  pybind11::handle get_comm() const {
    return pybind11::handle(PyMPIComm_New(ptr->get_comm()));
  }
  
  int get_nprocs() const { return ptr->get_nprocs(); }
  int get_nprow() const { return ptr->get_nprow(); }
  int get_npcol() const { return ptr->get_npcol(); }

  std::tuple<int,int> get_shape() const {
    return ptr->get_size();
  }
  
  int get_myrank() const { return ptr->get_myrank(); }
  int get_myrow() const { return ptr->get_myrow(); }
  int get_mycol() const { return ptr->get_mycol(); }

  std::tuple<int,int> get_mine() const {
    return ptr->get_my_coordinate();
  }
  
  bool is_row_major() const { return ptr->is_row_major(); }

  bool is_col_major() const { return ptr->is_col_major(); }

private:
  grid* ptr;
};

} // end namespace rokko

#endif // PYROKKO_GRID_HPP
