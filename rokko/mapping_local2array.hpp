/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MAPPING_LOCAL2ARRAY_HPP
#define ROKKO_MAPPING_LOCAL2ARRAY_HPP

#include <rokko/matrix_major.hpp>
#include <rokko/mapping_common_sizes.hpp>

namespace rokko {

template <typename MATRIX_MAJOR>
class mapping_local2array : virtual public mapping_common_sizes {
public:
  explicit mapping_local2array() {
    set_default_lld();
    set_default_length_array();
  }
  explicit mapping_local2array(int lld_in) : lld(lld_in) {
    set_default_length_array();
  }
  explicit mapping_local2array(int lld_in, int lld2) : lld(lld_in) {
    set_length_array(lld * lld2);
  }

  void set_length_array(int value) { length_array = value; }
  int get_length_array() const { return length_array; }

  int get_lld() const { return lld; };

  void set_lld(int value) { lld = value; };
  void set_default_lld() { set_lld(get_default_lld()); }

  int get_default_lld() const;

  int get_default_length_array() const;

  void set_default_length_array() { set_length_array(get_default_length_array()); }

  int get_array_index(int local_i, int local_j) const;

  bool is_row_major() const {
    return std::is_same<MATRIX_MAJOR, matrix_row_major>::value;
  }
  bool is_col_major() const {
    return std::is_same<MATRIX_MAJOR, matrix_col_major>::value;
  }

private:
  int length_array;
  int lld;
};


template<>
inline int mapping_local2array<rokko::matrix_row_major>::get_default_lld() const {
  return get_n_local();
}

template<>
inline int mapping_local2array<rokko::matrix_col_major>::get_default_lld() const {
  return get_m_local();
}

template<>
inline int mapping_local2array<rokko::matrix_row_major>::get_default_length_array() const {
  return get_m_local() * lld;
}

template<>
inline int mapping_local2array<rokko::matrix_col_major>::get_default_length_array() const {
  return lld * get_n_local();
}

template<>
inline int mapping_local2array<rokko::matrix_row_major>::get_array_index(int local_i, int local_j) const {
  return local_i * lld + local_j;
}

template<>
inline int mapping_local2array<rokko::matrix_col_major>::get_array_index(int local_i, int local_j) const {
  return  local_i + local_j * lld;
}

} // namespace rokko

#endif // ROKKO_MAPPING_LOCAL2ARRAY_HPP
