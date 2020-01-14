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

#ifndef PYROKKO_DISTRIBUTED_CRS_MATRIX_HPP
#define PYROKKO_DISTRIBUTED_CRS_MATRIX_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <rokko/pyrokko_communicator.hpp>
#include <rokko/pyrokko_mapping_1d.hpp>
#include <rokko/distributed_crs_matrix.hpp>
#include <rokko/utility/tuple_to_array.hpp>

namespace rokko {

class wrap_distributed_crs_matrix : public rokko::distributed_crs_matrix {
public:
  explicit wrap_distributed_crs_matrix(wrap_mapping_1d const& map, int num_entries_per_row)
    : distributed_crs_matrix(static_cast<mapping_1d>(map), num_entries_per_row) {}

  template<typename SOLVER>
  explicit wrap_distributed_crs_matrix(std::tuple<int,int> const& dims, SOLVER& solver_in)
    : distributed_crs_matrix(rokko::to_array(dims), solver_in) {}

  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) const {
    rokko::distributed_crs_matrix::insert(row, cols, values);
  }

  ~wrap_distributed_crs_matrix() {}

};


} // end namespace rokko

#endif // PYROKKO_DISTRIBUTED_CRS_MATRIX_HPP
