/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MATRIX_COMMON_HPP
#define ROKKO_MATRIX_COMMON_HPP

namespace rokko {

class matrix_common_sizes {
public:
  int get_m_local() const { return m_local; }
  int get_n_local() const { return n_local; }
  void set_m_local(int m_local_in) { m_local = m_local_in; }
  void set_n_local(int n_local_in) { n_local = n_local_in; }
  int get_nprow() const { return nprow; }
  int get_npcol() const { return npcol; }
  int get_nprocs() const { return nprocs; }
  int get_myrank() const { return myrank; }
  int get_myrow() const { return myrow; }
  int get_mycol() const { return mycol; }
  grid const& get_grid() const { return g; }

protected:
  explicit matrix_common_sizes() {}
  /*explicit matrix_common_sizes(matrix_common_sizes const& a)
    : m_local(a.m_local), n_local(a.n_local), myrank(a.myrank) {}*/
  explicit matrix_common_sizes(grid const& g_in) : g(g_in) {}
  //int get_dim() const { return dim_; }
  //int get_dim_local() const { return dim_local; }
  //int get_block_size() const { return block_size; }

  int m_local, n_local;
  grid g;
  // common variables of class grid
  int myrank, nprocs;
  int myrow, mycol;
  int nprow, npcol;
private:

};

} // namespace rokko

#endif // ROKKO_MATRIX_COMMON_HPP
