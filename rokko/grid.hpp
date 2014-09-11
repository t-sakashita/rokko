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
  // This constructor produces 2D grid whose row size is longer.  In
  // EigenK's grid creation, the variables nprow, npcol are changed
  // with this constructor's one, and in EigenK's sample main program,
  // some strange substitution is done: "NPROW = size_of_col NPCOL =
  // size_of_row".  So, this constructor's 2D grid size is consistent
  // with EigenK's 2D grid one.
  explicit grid(MPI_Comm comm_in = MPI_COMM_WORLD) : comm(comm_in) {
    initialize(grid_row_major);
  }
  template <typename GRID_MAJOR>
  grid(MPI_Comm comm_in, GRID_MAJOR const& grid_major) : comm(comm_in) {
    initialize(grid_major);
  }
  // The following 2D grid size is consistent with EigenK's 2D grid
  // creator's one, but it is actually changed with actual EigenK's 2D
  // grid row/col size. So we don't use it.
  /*
  grid(MPI_Comm comm_in = MPI_COMM_WORLD) {
    comm = comm_in;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myrank);

    npcol = int(std::sqrt(nprocs + 0.5));
    while (1) {
      if ( npcol == 1 ) break;
      if ( (nprocs % npcol) == 0 ) break;
      npcol = npcol - 1;
    }
    nprow = nprocs / npcol;
    myrow = calculate_grid_row(myrank);
    mycol = calculate_grid_col(myrank);
  }
  */
  MPI_Comm get_comm() const { return comm; }
  int get_nprocs() const { return nprocs; }
  int get_nprow() const { return nprow; }
  int get_npcol() const { return npcol; }
  int get_myrank() const { return myrank; }
  int get_myrow() const { return myrow; }
  int get_mycol() const { return mycol; }

  bool is_row_major() const { return is_row; }
  bool is_col_major() const { return !is_row_major(); }

  int calculate_grid_row(int proc_rank) const { 
    return is_row ? proc_rank / npcol : proc_rank % nprow;
  }
  int calculate_grid_col(int proc_rank) const { 
    return is_row ? proc_rank % npcol : proc_rank / nprow;
  }
protected:
  template <typename GRID_MAJOR>
  void initialize(GRID_MAJOR) {
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
