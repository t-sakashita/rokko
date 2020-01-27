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

#ifndef PYROKKO_DISTRIBUTED_MFREE_HPP
#define PYROKKO_DISTRIBUTED_MFREE_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <rokko/distributed_mfree.hpp>
#include <rokko/eigen3.hpp>

#include <rokko/mpi_communicator.hpp>
#include <rokko/pyrokko_communicator.hpp>

namespace rokko {

class wrap_distributed_mfree : public distributed_mfree_holder {
public:
  wrap_distributed_mfree(std::function<void(ConstMapVec,MapVec)> const& multiply, int dim)
    : distributed_mfree_holder(multiply, skel::mapping_1d(dim, rokko::mpi_comm{MPI_COMM_WORLD})) {}

  wrap_distributed_mfree(std::function<void(ConstMapVec,MapVec)> const& multiply, int dim, pybind11::handle const& comm_handle)
    : distributed_mfree_holder(multiply, skel::mapping_1d(dim, wrap_communicator{comm_handle})) {}

  ~wrap_distributed_mfree() = default;
};


} // end namespace rokko

#endif // PYROKKO_DISTRIBUTED_MFREE_HPP
