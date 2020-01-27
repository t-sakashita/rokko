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

#ifndef ROKKO_MAPPING_1D_COMMON_HPP
#define ROKKO_MAPPING_1D_COMMON_HPP

namespace rokko {

namespace detail {

class mapping_1d_common {
public:
  mapping_1d_common() = default;
  virtual ~mapping_1d_common() = default;

  virtual int get_dim() const = 0;
  virtual int get_num_local_rows() const = 0;
};

} // end namespace detail

} // end namespace rokko

#endif // ROKKO_MAPPING_1D_COMMON_HPP
