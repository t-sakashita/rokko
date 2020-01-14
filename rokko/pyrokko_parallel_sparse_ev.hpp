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

#ifndef PYROKKO_PARALLEL_SPARSE_EV_HPP
#define PYROKKO_PARALLEL_SPARSE_EV_HPP

#include <rokko/parallel_sparse_ev.hpp>
#include <rokko/pyrokko_mapping_1d.hpp>
#include <rokko/pyrokko_distributed_crs_matrix.hpp>
#include <rokko/pyrokko_distributed_mfree.hpp>

namespace rokko {

class wrap_parallel_sparse_ev : public parallel_sparse_ev {
public:
  wrap_parallel_sparse_ev(std::string const& solver_name) : parallel_sparse_ev(solver_name) {}
  
  wrap_parallel_sparse_ev() {}

  void initialize() {
    int num = 1;
    char** ptr = NULL;
    parallel_sparse_ev::initialize(num, ptr);
  }

  wrap_mapping_1d default_mapping(int dim, pybind11::handle const& comm_handle) const {
    return static_cast<wrap_mapping_1d>(parallel_sparse_ev::default_mapping(dim, static_cast<mpi_comm>(wrap_communicator{comm_handle})));
  }

  wrap_parameters diagonalize(wrap_distributed_crs_matrix& mat, wrap_parameters const& params) {
    return parallel_sparse_ev::diagonalize(mat, parameters(params));
  }

  wrap_parameters diagonalize(distributed_mfree* mat, wrap_parameters const& params) {
    return parallel_sparse_ev::diagonalize(*mat, parameters(params));
  }

  wrap_parameters diagonalize(wrap_distributed_mfree& mat, wrap_parameters const& params) {
    return parallel_sparse_ev::diagonalize(mat, parameters(params));
  }

  std::vector<double> python_eigenvector(int k) const {
    std::vector<double> vec;
    parallel_sparse_ev::eigenvector(k, vec);
    return vec;
  }
};

} // end namespace rokko

#endif // PYROKKO_PARALLEL_SPARSE_EV_HPP
