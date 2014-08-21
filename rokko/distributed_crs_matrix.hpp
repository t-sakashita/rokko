/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_DISTRIBUTED_CRS_MATRIX_HPP
#define ROKKO_DISTRIBUTED_CRS_MATRIX_HPP

#include <rokko/mapping_1d.hpp>

namespace rokko {

namespace detail {

class distributed_crs_matrix_base {
public:
  virtual void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) = 0;
  virtual void complete() = 0;
};

} // end namespace detail

class distributed_crs_matrix {
public:
  template<typename SOLVER>
  distributed_crs_matrix(mapping_1d const& map, SOLVER const& solver_in) {
    mat = solver_in.create_distributed_crs_matrix(map);
  }
  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) {
    mat->insert(row, cols, values);
  }
  void complete() {
    mat->complete();
  }
  detail::distributed_crs_matrix_base* get_matrix() {
    return mat;
  }
  //private:
  detail::distributed_crs_matrix_base* mat;
};

} // end namespace rokko

#endif // ROKKO_DISTRIBUTED_CRS_MATRIX_HPP

