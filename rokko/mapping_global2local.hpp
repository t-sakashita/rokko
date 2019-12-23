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

#ifndef ROKKO_MAPPING_GLOBAL2LOCAL_HPP
#define ROKKO_MAPPING_GLOBAL2LOCAL_HPP

#include <array>
#include <rokko/grid.hpp>
#include <rokko/mapping_common_sizes.hpp>

namespace rokko {

class mapping_global2local : virtual public mapping_common_sizes {
public:
  explicit mapping_global2local() {}

  explicit mapping_global2local(std::array<int,2> global_size_in, std::array<int,2> block_size_in, grid const& g_in)
    : global_size(global_size_in), block_size(block_size_in),
      g(g_in), my_coordinate(g_in.get_my_coordinate()), grid_size(g_in.get_size()) {
    set_default_local_size();
    set_stride();
  }

  explicit mapping_global2local(int global_dim, int block_size, grid const& g_in)
    : mapping_global2local({global_dim, global_dim}, {block_size, block_size}, g_in) {}

  std::array<int,2> get_stride_mine() const { return stride_mine; }
  std::array<int,2> get_stride_grid() const { return stride_grid; }

  void set_stride_mine() {
    stride_mine[0] = my_coordinate[0] * block_size[0];
    stride_mine[1] = my_coordinate[1] * block_size[1];
  }

  void set_stride_grid() {
    stride_grid[0] = block_size[0] * (grid_size[0] - 1);
    stride_grid[1] = block_size[1] * (grid_size[1] - 1);
  }

  void set_stride() {
    set_stride_mine();
    set_stride_grid();
  }

  void set_block_size(std::array<int,2> block_size_in) {
    block_size = block_size_in;
  }

  int get_mb() const { return std::get<0>(block_size); }
  int get_nb() const { return std::get<1>(block_size); }
  std::array<int,2> get_block_size() const { return block_size; }

  int get_m_global() const { return std::get<0>(global_size); }
  int get_n_global() const { return std::get<1>(global_size); }
  std::array<int,2> get_global_size() const { return global_size; }

  void set_default_local_size() {
    set_local_size(calculate_default_local_size(my_coordinate));
  }

  template <int IND>
  int calculate_default_local_size(int proc) const {
    const int quotient = global_size[IND] / block_size[IND];
    const int local_num_block = (quotient - 1 - proc) / grid_size[IND] + 1;
    const int remainder_block_proc = quotient % grid_size[IND];
    const int local_remainder_block = (proc == remainder_block_proc) ? global_size[IND] % block_size[IND] : 0;  // size of a remainder block (< block_size)
    return local_num_block * block_size[IND] + local_remainder_block;
  }

  std::array<int,2> calculate_default_local_size(std::array<int,2> proc) const {
    return {calculate_default_local_size<0>(proc[0]), calculate_default_local_size<1>(proc[1])};
  }

  template <int IND>
  int translate_l2g(int local_ind) const {
    return stride_mine[IND] + local_ind + (local_ind / block_size[IND]) * stride_grid[IND];
  }

  std::array<int,2> translate_l2g(std::array<int,2> local_indices) const {
    return {translate_l2g<0>(local_indices[0]), translate_l2g<1>(local_indices[1])};
  }

  int translate_l2g_row(int local_i) const {
    return translate_l2g<0>(local_i);
  }

  int translate_l2g_col(int local_j) const {
    return translate_l2g<1>(local_j);
  }

  template <int IND>
  int translate_g2l(int global_index) const {
    const int local_offset_block = global_index / block_size[IND];
    return (local_offset_block - my_coordinate[IND]) / grid_size[IND] * block_size[IND] + global_index % block_size[IND];
  }

  std::array<int,2> translate_g2l(std::array<int,2> global_indices) const {
    return {translate_g2l<0>(global_indices[0]), translate_g2l<1>(global_indices[1])};
  }

  int translate_g2l_row(int global_i) const {
    return translate_g2l<0>(global_i);
  }

  int translate_g2l_col(int global_j) const {
    return translate_g2l<1>(global_j);
  }

  template <int IND>
  bool is_gindex(int global_index) const {
    const int local_offset_block = global_index / block_size[IND];
    return (local_offset_block % grid_size[IND]) == my_coordinate[IND];
  }

  bool is_gindex_myrow(int global_i) const {
    return is_gindex<0>(global_i);
  }

  bool is_gindex_mycol(int global_j) const {
    return is_gindex<1>(global_j);
  }

  bool is_gindex(std::array<int,2> global_indices) const {
    return is_gindex<0>(global_indices[0]) && is_gindex<1>(global_indices[1]);
  }

  bool is_gindex(int global_i, int global_j) const {
    return is_gindex<0>(global_i) && is_gindex_mycol(global_j);
  }

  MPI_Comm get_comm() const { return g.get_comm(); }
  int get_nprow() const { return g.get_nprow(); }
  int get_npcol() const { return g.get_npcol(); }
  int get_nprocs() const { return g.get_nprocs(); }
  int get_myrank() const { return g.get_myrank(); }
  int get_myrow() const { return g.get_myrow(); }
  int get_mycol() const { return g.get_mycol(); }
  grid const& get_grid() const { return g; }

private:
  std::array<int,2> global_size;
  std::array<int,2> block_size;
  std::array<int,2> stride_mine, stride_grid;
  grid g;
  // common variables of class grid
  std::array<int,2> grid_size;
  std::array<int,2> my_coordinate;

  // block_number is also needed?
};

} // namespace rokko

#endif // ROKKO_MAPPING_GLOBAL2LOCAL_HPP
