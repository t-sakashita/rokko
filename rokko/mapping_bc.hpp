/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2014 by Synge Todo <wistaria@comp-phys.org>,
*                            Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MAPPING_BC_HPP
#define ROKKO_MAPPING_BC_HPP

#include <rokko/grid.hpp>

namespace rokko {

class mapping_bc {
public:
  explicit mapping_bc() : dim_(0), g_() {}
  explicit mapping_bc(grid const& g, int dim, int block_size) :
    g_(g), dim_(dim), block_size_(block_size) {
    dim_local_ = dim / block_size;
  }
  int get_dim() const { return dim_; }
  int get_dim_local() const { return dim_local_; }
  int get_block_size() const { return block_size_; }
  int get_lld() const { return lld; }
  void set_lld(int value) { lld = value; }
  grid const& get_grid() const { return g_; }

private:
  int dim_;
  int block_size_;
  // block_number is also needed?
  int dim_local_;
  int lld;
  grid g_;
};

} // namespace rokko

#endif // ROKKO_MAPPING_BC_HPP
