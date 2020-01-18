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


namespace rokko {

using MapVec = Eigen::Map<Eigen::Vector<double>>;
using ConstMapVec = const Eigen::Map<const Eigen::Vector<double>>;

class wrap_distributed_mfree : public rokko::distributed_mfree_default {
public:
  wrap_distributed_mfree(std::function<void(ConstMapVec,MapVec)> const& multiply, int dim, MPI_Comm comm = MPI_COMM_WORLD)
    : multiply_(multiply), rokko::distributed_mfree_default(dim, rokko::mpi_comm{comm}) {}
  ~wrap_distributed_mfree() {}

  void multiply(const double* x, double* y) const {
    int num_local_rows = get_num_local_rows();
    ConstMapVec  X(x, num_local_rows);
    MapVec  Y(y, num_local_rows);
    multiply_(X, Y);
  }

private:
  std::function<void(ConstMapVec,MapVec)> multiply_;
};


} // end namespace rokko

#endif // PYROKKO_DISTRIBUTED_MFREE_HPP
