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

#ifndef ROKKO_DISTRIBUTED_MFREE_HPP
#define ROKKO_DISTRIBUTED_MFREE_HPP

#include <rokko/mapping_1d.hpp>

namespace rokko {

class distributed_mfree {
public:
  distributed_mfree() {}
  ~distributed_mfree() {}
  mapping_1d get_mapping_1d() {
    return map_;
  }
  virtual void multiply(const double* x, double* y) const = 0;
protected:
  mapping_1d map_;
};

} // end namespace rokko

#endif // ROKKO_DISTRIBUTED_MFREE_HPP
