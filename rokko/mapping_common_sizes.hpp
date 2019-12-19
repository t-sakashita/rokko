/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MAPPING_COMMON_SIZES_HPP
#define ROKKO_MAPPING_COMMON_SIZES_HPP

namespace rokko {

class mapping_common_sizes {
public:
  int get_m_local() const { return std::get<0>(local_size); }
  int get_n_local() const { return std::get<1>(local_size); }
  void set_m_local(int m_local) { std::get<0>(local_size) = m_local; }
  void set_n_local(int n_local) { std::get<1>(local_size) = n_local; }

protected:
  explicit mapping_common_sizes() {}

private:
  std::array<int,2> local_size;
};

} // namespace rokko

#endif // ROKKO_MAPPING_COMMON_SIZES_HPP
