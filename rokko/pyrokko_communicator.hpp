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

#ifndef PYROKKO_COMMUNICATOR_HPP
#define PYROKKO_COMMUNICATOR_HPP

#include <mpi4py/mpi4py.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <rokko/mpi_communicator.hpp>

namespace rokko {

namespace py = pybind11;

class wrap_communicator : public mpi_comm {
public:
  using mpi_comm::mpi_comm;
  wrap_communicator(pybind11::handle const& comm_handle) {
    const int rc = import_mpi4py();
    assert(rc==0);
    set_comm(*PyMPIComm_Get(comm_handle.ptr()));
  }
};

} // end namespace rokko

#endif // PYROKKO_COMMUNICATOR_HPP
