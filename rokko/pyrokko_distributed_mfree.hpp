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
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include <rokko/eigen3.hpp>

#include <rokko/mpi_communicator.hpp>
#include <rokko/distributed_mfree.hpp>
#include <rokko/pyrokko_communicator.hpp>
#include <rokko/pyrokko_mapping_1d.hpp>

namespace rokko {

class wrap_distributed_mfree : public distributed_mfree_holder {
public:
  wrap_distributed_mfree(std::function<void(ConstMapVec,MapVec)> const& multiply, int dim)
    : distributed_mfree_holder(multiply, skel::mapping_1d(dim, rokko::mpi_comm{MPI_COMM_WORLD})) {}

  wrap_distributed_mfree(std::function<void(ConstMapVec,MapVec)> const& multiply, int dim, pybind11::handle const& comm_handle)
    : distributed_mfree_holder(multiply, skel::mapping_1d(dim, wrap_communicator{comm_handle})) {}

  wrap_distributed_mfree(std::function<void(ConstMapVec,MapVec)> const& multiply, wrap_mapping_1d const& map)
    : distributed_mfree_holder(multiply, skel::mapping_1d(map.get_dim(), map.get_mpi_comm())) {}

  ~wrap_distributed_mfree() = default;
};


class distributed_mfree_inherit : public distributed_mfree_default {
public:
  using distributed_mfree_default::distributed_mfree_default;

  distributed_mfree_inherit(int dim, pybind11::handle const& comm_handle)
  : distributed_mfree_default(skel::mapping_1d(dim, wrap_communicator{comm_handle})) {}

  virtual void multiply(const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> y) const { // Not pure virtual to display error message
    std::cerr << "distributed_mfree_inherit : multiply is not implemented in Python" << std::endl;
  };

  void multiply(const double* x, double *const y) const override {
    const auto num_local_rows = get_num_local_rows();
    Eigen::Vector<double> X, Y;
    new (&X) Eigen::Map<const Eigen::Vector<double>>(x, num_local_rows);
    new (&Y) Eigen::Map<Eigen::Vector<double>>(y, num_local_rows);
    multiply(X, Y);
    new (&X) Eigen::Vector<double>{};
    new (&Y) Eigen::Vector<double>{};
  }
};

class wrap_distributed_mfree_inherit : public distributed_mfree_inherit {
public:
  using distributed_mfree_inherit::distributed_mfree_inherit;

  void multiply(const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> y) const override {
    PYBIND11_OVERLOAD(void, distributed_mfree_inherit, multiply, x, y);
  }
};

} // end namespace rokko
