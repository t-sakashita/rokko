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
      g(g_in), myrank(g_in.get_myrank()), nprocs(g_in.get_nprocs()),
      myrow(g_in.get_myrow()), mycol(g_in.get_mycol()),
      nprow(g_in.get_nprow()), npcol(g_in.get_npcol()) {
    set_default_local_size();
    set_stride();
  }

  explicit mapping_global2local(int global_dim, int block_size, grid const& g_in)
    : mapping_global2local({global_dim, global_dim}, {block_size, block_size}, g_in) {}

  /*explicit mapping_global2local(int m_global_in, int n_global_in, int m_local_in, int n_local_in, grid const& g_in) : mapping_common_sizes(g_in), m_global(m_global_in), n_global(n_global_in), m_local(m_local_in), n_local(n_local_in) {
  // get grid information
  myrank = g.get_myrank(); nprocs = g.get_nprocs();
  myrow = g.get_myrow(); mycol = g.get_mycol();
  nprow = g.get_nprow(); npcol = g.get_npcol();
  mb = m_global / m_local;
  nb = n_global / n_local;
  set_stride();
  }*/

  int get_stride_myrow() const { return stride_mine[0]; }
  int get_stride_mycol() const { return stride_mine[1]; }

  void set_stride() {
    stride_mine[0] = myrow * get_mb();
    stride_grid[0] = get_mb() * (nprow - 1);
    stride_mine[1] = mycol * get_nb();
    stride_grid[1] = get_nb() * (npcol - 1);
  }
  void set_block_size(std::array<int,2> block_size_in) {
    block_size = block_size_in;
  }

  int get_mb() const { return std::get<0>(block_size); }
  int get_nb() const { return std::get<1>(block_size); }

  int get_m_global() const { return std::get<0>(global_size); }
  int get_n_global() const { return std::get<1>(global_size); }

  void set_default_local_size() {
    set_local_size({calculate_row_size(), calculate_col_size()});
  }

  int calculate_row_size(int proc_row) const {
    int tmp = get_m_global() / get_mb();
    int local_num_block_rows = (tmp - proc_row -1) / nprow + 1;
    int rest_block_row = tmp % nprow; // size of a residue block (< mb)
    int local_rest_block_rows = (proc_row == rest_block_row) ? get_m_global() % get_mb() : 0;

    return  local_num_block_rows * get_mb() + local_rest_block_rows;
  }

  int calculate_row_size() const {
    return calculate_row_size(myrow);
  }

  int calculate_col_size(int proc_col) const {
    int tmp = get_n_global() / get_nb();
    int local_num_block_cols = (tmp - proc_col -1) / npcol + 1;
    int rest_block_col = tmp % npcol; // size of a residue block (< nb)
    int local_rest_block_cols = (proc_col == rest_block_col) ? get_n_global() % get_nb() : 0;
    return local_num_block_cols * get_nb() + local_rest_block_cols;
  }

  int calculate_col_size() const {
    return calculate_col_size(mycol);
  }

  int translate_l2g_row(const int& local_i) const {
    return stride_mine[0] + local_i + (local_i / get_mb()) * stride_grid[0];
  }

  int translate_l2g_col(const int& local_j) const {
    return stride_mine[1] + local_j + (local_j / get_nb()) * stride_grid[1];
  }

  int translate_g2l_row(const int& global_i) const {
    int local_offset_block = global_i / get_mb();
    return (local_offset_block - myrow) / nprow * get_mb() + global_i % get_mb();
  }

  int translate_g2l_col(const int& global_j) const {
    const int local_offset_block = global_j / get_nb();
    return (local_offset_block - mycol) / npcol * get_nb() + global_j % get_nb();
  }

  bool is_gindex_myrow(const int& global_i) const {
    int local_offset_block = global_i / get_mb();
    return (local_offset_block % nprow) == myrow;
  }

  bool is_gindex_mycol(const int& global_j) const {
    int local_offset_block = global_j / get_nb();
    return (local_offset_block % npcol) == mycol;
  }

  bool is_gindex(const int& global_i, const int& global_j) const {
    return is_gindex_myrow(global_i) && is_gindex_mycol(global_j);
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
  int myrank, nprocs;
  int myrow, mycol;
  int nprow, npcol;
  // block_number is also needed?
};

} // namespace rokko

#endif // ROKKO_MAPPING_GLOBAL2LOCAL_HPP
