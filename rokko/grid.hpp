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

#ifndef ROKKO_GRID_HPP
#define ROKKO_GRID_HPP

#include <mpi.h>
#include <cmath>
#include <boost/type_traits/is_same.hpp>

namespace rokko {

extern struct grid_row_major_t {} grid_row_major;

extern struct grid_col_major_t {} grid_col_major;

class grid {
public:
  /*grid(grid const& g_in)
    : comm(g_in.get_comm()), nprocs(g_in.get_nprocs()), nprow(g_in.get_nprow()), npcol(g_in.get_npcol()), myrank(g_in.get_myrank()), myrow(g_in.get_myrow()), mycol(g_in.get_mycol()), is_row(g_in.is_row_major()) {}*/
  explicit grid(MPI_Comm comm_in = MPI_COMM_WORLD) : comm(comm_in) {
    initialize(grid_row_major);
  }
  template <typename GRID_MAJOR>
  grid(MPI_Comm& comm_in, GRID_MAJOR const& grid_major = grid_row_major) : comm(comm_in) {
    initialize(grid_major);
  }

  MPI_Comm get_comm() const { return comm; }
  int get_nprocs() const { return nprocs; }
  int get_nprow() const { return nprow; }
  int get_npcol() const { return npcol; }
  int get_myrank() const { return myrank; }
  int get_myrow() const { return myrow; }
  int get_mycol() const { return mycol; }

  bool is_row_major() const { return is_row; }
  bool is_col_major() const { return !is_row; }

  int calculate_grid_row(int proc_rank) const { 
    return is_row ? proc_rank / npcol : proc_rank % nprow;
  }
  int calculate_grid_col(int proc_rank) const { 
    return is_row ? proc_rank % npcol : proc_rank / nprow;
  }
protected:
  template <typename GRID_MAJOR>
  void initialize(GRID_MAJOR) {
    if (comm == MPI_COMM_NULL) {
      myrank = -1;
      myrow = -1;
      mycol = -1;
      nprocs = 0;
      nprow = 0; // to avoid zero dividing and blacs does not take (0,0) size grid
      npcol = 0; // to avoid zero dividing and blacs does not take (0,0) size grid
      return;
    }
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myrank);
    is_row = boost::is_same<GRID_MAJOR, grid_row_major_t>::value;
    nprow = int(std::sqrt(nprocs + 0.5));
    while (1) {
      if ( nprow == 1 ) break;
      if ( (nprocs % nprow) == 0 ) break;
      nprow = nprow - 1;
    }
    npcol = nprocs / nprow;
    myrow = calculate_grid_row(myrank);
    mycol = calculate_grid_col(myrank);
  }
private:
  MPI_Comm comm;
  int nprocs, myrank;
  int nprow, npcol;
  int myrow, mycol;
  bool is_row;
};

} // namespace rokko

#endif // ROKKO_GRID_HPP
