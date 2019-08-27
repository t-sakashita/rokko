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

#ifndef PYROKKO_COMMUNICATOR_HPP
#define PYROKKO_COMMUNICATOR_HPP

#include <mpi4py/mpi4py.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace rokko {

namespace py = pybind11;

class wrap_communicator {
public:
  wrap_communicator() : comm(MPI_COMM_WORLD) {}
  
  wrap_communicator(MPI_Comm const& comm) : comm(comm) {}
  
  MPI_Comm get_comm() const { return comm; }

private:
  MPI_Comm comm;
};

} // end namespace rokko

#endif // PYROKKO_COMMUNICATOR_HPP
