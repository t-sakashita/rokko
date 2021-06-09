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

#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <rokko/mapping_bc.hpp>
#include <rokko/pyrokko_grid.hpp>
#include <rokko/pyrokko_matrix_major_enum.hpp>
#include <rokko/utility/tuple_to_array.hpp>

namespace rokko {

class base_mapping_bc {
public:
  virtual ~base_mapping_bc() = default;
};

template<typename MATRIX_MAJOR>
class wrap_mapping_bc : public mapping_bc<MATRIX_MAJOR>, public base_mapping_bc {
public:
  wrap_mapping_bc(std::tuple<int,int> const& global_size, std::tuple<int,int> const& block_size, wrap_grid const& g)
    : mapping_bc<MATRIX_MAJOR>(to_array(global_size), to_array(block_size), g) {}

  wrap_mapping_bc(int global_dim, int block_size, int lld, wrap_grid const& g)
    : mapping_bc<MATRIX_MAJOR>(global_dim, block_size, lld, g) {}

  wrap_mapping_bc(int global_dim, int block_size, wrap_grid const& g)
    : mapping_bc<MATRIX_MAJOR>(global_dim, block_size, g) {}

  wrap_mapping_bc(mapping_bc<MATRIX_MAJOR> const& map) : mapping_common_sizes(map), mapping_bc<MATRIX_MAJOR>(map) {}

  bool has_global_indices(std::tuple<int,int> const& global) const {
    return mapping_bc<MATRIX_MAJOR>::has_global_indices(to_array(global));
  }

  auto get_mb() const {
    return mapping_bc<MATRIX_MAJOR>::get_mb();
  }

  auto get_nb() const {
    return mapping_bc<MATRIX_MAJOR>::get_nb();
  }

  auto get_m_local() const {
    return mapping_bc<MATRIX_MAJOR>::get_m_local();
  }

  auto get_n_local() const {
    return mapping_bc<MATRIX_MAJOR>::get_n_local();
  }

  std::tuple<int,int> get_block_shape() const {
    return mapping_bc<MATRIX_MAJOR>::get_block_size();
  }
  
  std::tuple<int,int> get_global_shape() const {
    return mapping_bc<MATRIX_MAJOR>::get_global_size();
  }
  
  std::tuple<int,int> get_local_shape() const {
    return mapping_bc<MATRIX_MAJOR>::get_local_size();
  }

  std::tuple<int,int> translate_l2g(std::tuple<int,int> const& local) const {
    return mapping_bc<MATRIX_MAJOR>::translate_l2g(to_array(local));
  }

  std::tuple<int,int> translate_g2l(std::tuple<int,int> const& global) const {
    return mapping_bc<MATRIX_MAJOR>::translate_g2l(to_array(global));
  }

  wrap_grid get_grid() const {
    return wrap_grid(get_grid());
  }

  std::string get_major_string() const {
    return mapping_bc<MATRIX_MAJOR>::is_col_major() ? "col" : "row";
  }
};

std::shared_ptr<base_mapping_bc> create_mapping_bc(std::tuple<int,int> const& global_size, std::tuple<int,int> const& block_size, wrap_grid const& wrap_g, matrix_major_enum const& major = col) {
  if (major == col)
    return std::make_shared<wrap_mapping_bc<matrix_col_major>>(global_size, block_size, wrap_g);
  else
    return std::make_shared<wrap_mapping_bc<matrix_row_major>>(global_size, block_size, wrap_g);
}

std::shared_ptr<base_mapping_bc> create_mapping_bc(int global_size, int block_size, int lld, wrap_grid const& wrap_g, matrix_major_enum const& major = col) {
  if (major == col)
    return std::make_shared<wrap_mapping_bc<matrix_col_major>>(global_size, block_size, lld, wrap_g);
  else
    return std::make_shared<wrap_mapping_bc<matrix_row_major>>(global_size, block_size, lld, wrap_g);
}

std::shared_ptr<base_mapping_bc> create_mapping_bc(int global_size, int block_size, wrap_grid const& wrap_g, matrix_major_enum const& major = col) {
  if (major == col)
    return std::make_shared<wrap_mapping_bc<matrix_col_major>>(global_size, block_size, wrap_g);
  else
    return std::make_shared<wrap_mapping_bc<matrix_row_major>>(global_size, block_size, wrap_g);
}

} // end namespace rokko
