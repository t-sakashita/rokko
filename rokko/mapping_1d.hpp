/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MAPPING_1D_HPP
#define ROKKO_MAPPING_1D_HPP

#include <rokko/grid_1d.hpp>

namespace rokko {

namespace detail {

class mapping_1d_base {
public:
  virtual int get_dim() = 0;
  virtual int get_num_local_rows() = 0;
};

} // end namespace detail

class mapping_1d {
public:
  explicit mapping_1d() {}

  template<typename SOLVER>
  explicit mapping_1d(SOLVER& solver_in) {
    map = solver_in.create_mapping_1d();  
  }
  template<typename SOLVER>
  explicit mapping_1d(int dim, rokko::grid_1d const& g, SOLVER& solver_in) {
    map = solver_in.create_mapping_1d(dim, g);
  }
  int get_dim() {
    return map->get_dim();
  }
  int get_num_local_rows() const {
    return map->get_num_local_rows();
  }
  detail::mapping_1d_base* get_mapping_1d() {
    return map;
  }
//private:
  detail::mapping_1d_base* map;
};

} // namespace rokko

#endif // ROKKO_MAPPING_1D_HPP
