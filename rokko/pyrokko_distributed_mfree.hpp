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

#ifndef PYROKKO_DISTRIBUTED_MFREE_HPP
#define PYROKKO_DISTRIBUTED_MFREE_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <rokko/distributed_mfree.hpp>
#include <rokko/localized_vector.hpp>


namespace rokko {

using MapVec = Eigen::Map<Eigen::Vector<double>>;
using ConstMapVec = const Eigen::Map<const Eigen::Vector<double>>;

class wrap_distributed_mfree : public rokko::distributed_mfree {
public:
  wrap_distributed_mfree(std::function<void(ConstMapVec,MapVec)> const& multiply, int dim, int num_local_rows)
    : multiply_(multiply), dim_(dim), num_local_rows_(num_local_rows), local_offset_(0) {}
  ~wrap_distributed_mfree() {}

  void multiply(const double* x, double* y) const {
    ConstMapVec  X(x, num_local_rows_);
    MapVec  Y(y, num_local_rows_);
    multiply_(X, Y);
  }

  int get_dim() const { return dim_; }
  int get_local_offset() const { return local_offset_; }
  int get_num_local_rows() const { return num_local_rows_; }

private:
  std::function<void(ConstMapVec,MapVec)> multiply_;
  int dim_;
  int num_local_rows_;
  int local_offset_;
};


} // end namespace rokko

#endif // PYROKKO_DISTRIBUTED_MFREE_HPP
