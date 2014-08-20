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

#include <rokko/factory.hpp>
#include <rokko/mapping_1d.hpp>

namespace rokko {

class parallel_sparse_solver;
  
class distributed_crs_matrix {
public:
  template<typename SOLVER>
  distributed_crs_matrix(std::string const& solver_name, mapping_1d const& map, SOLVER const& solver_in) {
    matrix_impl_ = detail::dc_matrix_factory::instance()->make_product(solver_name);
    matrix_impl_->initialize(map);
  }
  distributed_crs_matrix() {
    matrix_impl_ = detail::dc_matrix_factory::instance()->make_product();
  }
  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) {
    matrix_impl_->insert(row, cols, values);
  }
  void complete() { matrix_impl_->complete(); }
// private:
  detail::dc_matrix_factory::product_pointer_type matrix_impl_;
};

} // end namespace rokko

#endif // ROKKO_DISTRIBUTED_CRS_MATRIX_HPP

