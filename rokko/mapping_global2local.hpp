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

#include <rokko/grid.hpp>
#include <rokko/matrix_common.hpp>

namespace rokko {

class mapping_global2local : virtual public matrix_common_sizes {
public:
  explicit mapping_global2local() {}
  explicit mapping_global2local(int m_global_in, int n_global_in, int block_size, grid const& g_in)
    : g(g_in), myrank(g_in.get_myrank()), nprocs(g_in.get_nprocs()), myrow(g_in.get_myrow()), mycol(g_in.get_mycol()), nprow(g_in.get_nprow()), npcol(g_in.get_npcol()),
      m_global(m_global_in), n_global(n_global_in), mb(block_size), nb(block_size) {
    set_default_local_size();
    set_stride();
  }

  /*explicit mapping_global2local(grid const& g_in, int m_global_in, int n_global_in, int m_local_in, int n_local_in) : matrix_common_sizes(g_in), m_global(m_global_in), n_global(n_global_in), m_local(m_local_in), n_local(n_local_in) {
    // get grid information
    myrank = g.get_myrank(); nprocs = g.get_nprocs();
    myrow = g.get_myrow(); mycol = g.get_mycol();
    nprow = g.get_nprow(); npcol = g.get_npcol();
    mb = m_global / m_local;
    nb = n_global / n_local;
    set_stride();
    }*/
    
  int get_stride_myrow() const { return stride_myrow; }
  int get_stride_mycol() const { return stride_mycol; }

  void set_stride() {
    stride_myrow = myrow * mb;
    stride_nprow = mb * (nprow - 1);
    stride_mycol = mycol * nb;
    stride_npcol = nb * (npcol - 1);    
  }
  void set_block_size(int mb_in, int nb_in) {
    mb = mb_in;
    nb = nb_in;
  }

  int get_mb() const { return mb; }
  int get_nb() const { return nb; }

  int get_m_global() const { return m_global; }
  int get_n_global() const { return n_global; }

  void set_local_size(int m_local_in, int n_local_in) {
    set_m_local(m_local_in);
    set_n_local(n_local = n_local_in);
  }

  void set_default_local_size() {
    set_local_size(calculate_row_size(), calculate_col_size());
  }

  int calculate_row_size(int proc_row) const {
    int tmp = m_global / mb;
    int local_num_block_rows = (tmp - proc_row -1) / nprow + 1;
    int rest_block_row = tmp % nprow; // size of a residue block (< mb)
    int local_rest_block_rows;
    if (proc_row == rest_block_row)
      local_rest_block_rows = m_global % mb;
    else
      local_rest_block_rows = 0;

    return  local_num_block_rows * mb + local_rest_block_rows;
  }

  int calculate_row_size() const {
    return calculate_row_size(myrow);
  }

  int calculate_col_size(int proc_col) const {
    int tmp = n_global / nb;
    int local_num_block_cols = (tmp - proc_col -1) / npcol + 1;
    int rest_block_col = tmp % npcol; // size of a residue block (< nb)
    int local_rest_block_cols;
    if (proc_col == rest_block_col) {
      local_rest_block_cols = n_global % nb;
    } else {
      local_rest_block_cols = 0;
    }
    return local_num_block_cols * nb + local_rest_block_cols;
  }

  int calculate_col_size() const {
    return calculate_col_size(mycol);
  }

  int translate_l2g_row(const int& local_i) const {
    return stride_myrow + local_i + (local_i / mb) * stride_nprow;
  }

  int translate_l2g_col(const int& local_j) const {
    return stride_mycol + local_j + (local_j / nb) * stride_npcol;
  }

  int translate_g2l_row(const int& global_i) const {
    int local_offset_block = global_i / mb;
    return (local_offset_block - myrow) / nprow * mb + global_i % mb;
  }

  int translate_g2l_col(const int& global_j) const {
    const int local_offset_block = global_j / nb;
    return (local_offset_block - mycol) / npcol * nb + global_j % nb;
  }

  bool is_gindex_myrow(const int& global_i) const {
    int local_offset_block = global_i / mb;
    return (local_offset_block % nprow) == myrow;
  }

  bool is_gindex_mycol(const int& global_j) const {
    int local_offset_block = global_j / nb;
    return (local_offset_block % npcol) == mycol;
  }

  bool is_gindex(const int& global_i, const int& global_j) const {
    return is_gindex_myrow(global_i) && is_gindex_mycol(global_j);
  }

  int get_nprow() const { return g.get_nprow(); }
  int get_npcol() const { return g.get_npcol(); }
  int get_nprocs() const { return g.get_nprocs(); }
  int get_myrank() const { return g.get_myrank(); }
  int get_myrow() const { return g.get_myrow(); }
  int get_mycol() const { return g.get_mycol(); }
  grid const& get_grid() const { return g; }

protected:
  int m_global, n_global;
  int mb, nb;
  int stride_myrow, stride_nprow, stride_mycol, stride_npcol;
  grid g;
  // common variables of class grid
  int myrank, nprocs;
  int myrow, mycol;
  int nprow, npcol;
private:
  // block_number is also needed?
};

} // namespace rokko

#endif // ROKKO_MAPPING_GLOBAL2LOCAL_HPP
