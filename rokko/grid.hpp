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

#ifndef ROKKO_GRID_HPP
#define ROKKO_GRID_HPP

#include <array>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <rokko/blacs.hpp>

namespace rokko {

extern struct grid_row_major_t {} grid_row_major;

extern struct grid_col_major_t {} grid_col_major;

class grid {
public:
  template <typename GRID_MAJOR>
  grid(MPI_Comm& comm_in, GRID_MAJOR const& grid_major = grid_row_major, int lld = 0)
    : comm(comm_in) {
    if (comm == MPI_COMM_NULL) {
      nprocs = 0;
      myrank = -1;
      set_size({0, 0});  // to avoid zero dividing and blacs does not take (0,0) size grid
      set_my_coordinate({-1,-1});
      return;
    }
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myrank);

    if (is_MPI_2dim_Cart(comm))
      calculate_sizes_cart(grid_major);
    else
      calculate_sizes_default(grid_major, lld);

    set_blacs_grid();
  }

  template <typename GRID_MAJOR>
  grid(MPI_Comm& comm_in, std::array<int,2> size_in, GRID_MAJOR const& grid_major = grid_row_major) : comm(comm_in) {
    MPI_Comm_size(comm, &nprocs);
    if ((std::get<0>(size_in) * std::get<1>(size_in)) != nprocs) {
      throw std::invalid_argument("grid::grid() : (rows * cols) != nprocs");
    }
    set_size(size_in);

    set_major<GRID_MAJOR>();

    set_my_coordinate();

    set_blacs_grid();
  }

  explicit grid(MPI_Comm comm_in = MPI_COMM_WORLD) : grid(comm_in, grid_row_major, 0) {}

  MPI_Comm get_comm() const { return comm; }
  int get_nprocs() const { return nprocs; }
  int get_nprow() const { return std::get<0>(size); }
  int get_npcol() const { return std::get<1>(size); }
  std::array<int,2> get_size() const { return size; }
  int get_myrank() const { return myrank; }
  int get_myrow() const { return std::get<0>(my_coordinate); }
  int get_mycol() const { return std::get<1>(my_coordinate); }
  std::array<int,2> get_my_coordinate() const { return my_coordinate; }
  int get_blacs_context() const { return blacs_context; }

  void set_size(std::array<int,2> size_in) { size = size_in; }
  void set_my_coordinate(std::array<int,2> my_coordinate_in) { my_coordinate = my_coordinate_in; }
  void set_my_coordinate() { set_my_coordinate(calculate_coordinate(myrank)); }

  template <typename GRID_MAJOR>
  void set_major() { is_row = std::is_same<GRID_MAJOR, grid_row_major_t>::value; }

  bool is_row_major() const { return is_row; }
  bool is_col_major() const { return !is_row; }

  int calculate_grid_row(int proc_rank) const { 
    return is_row ? proc_rank / get_npcol() : proc_rank % get_nprow();
  }
  int calculate_grid_col(int proc_rank) const { 
    return is_row ? proc_rank % get_npcol() : proc_rank / get_nprow();
  }

  std::array<int,2> calculate_coordinate(int proc_rank) const {
    return {calculate_grid_row(proc_rank), calculate_grid_col(proc_rank)};
  }

  int calculate_rank_form_coords(int proc_row, int proc_col) const {
    return is_row ? (proc_row * get_npcol() + proc_col)
      : (proc_col * get_nprow() + proc_row);
  }

protected:
  template <typename GRID_MAJOR>
  void calculate_sizes_default(GRID_MAJOR, int lld) {
    int nprow;
    if (lld > 0) {
      if ((nprocs % lld) != 0) {
        throw std::invalid_argument("The number of processes should be a multiple of lld.");
      }
      nprow = (is_row ? nprocs / lld : lld);
    } else {
      nprow = find_square_root_like_divisor(nprocs);
    }
    set_size({nprow, nprocs/nprow});
    set_major<GRID_MAJOR>();
    set_my_coordinate();
  }

  static int find_square_root_like_divisor(int n) {
    int i = int(std::sqrt((double)n));
    for (; i > 1; --i) {
      if ( (n % i) == 0 ) break;
    }
    return i;
  }

  void calculate_sizes_cart(grid_row_major_t) {
    constexpr int cart_dim = 2;
    std::array<int,2> dims, periods, coords;
    /* int ierr = */ MPI_Cart_get(comm, cart_dim, dims.data(), periods.data(), coords.data());
    set_size(dims);
    set_my_coordinate(coords);

    set_major<grid_row_major_t>();
  }

  static void calculate_sizes_cart(grid_col_major_t) {
    throw std::invalid_argument("calculate_sizes_cart : MPI Cartesian doesn't support grid_col_major.  Use it with grid_row_major.");
  }

  static bool is_MPI_Cart(MPI_Comm comm) {
    int topo_type;
    /* int ierr = */ MPI_Topo_test(comm, &topo_type);
    return topo_type == MPI_CART;
  }

  static bool is_MPI_2dim_Cart(MPI_Comm comm) {
    if (is_MPI_Cart(comm)) {
      int cart_dim;
      /* int ierr = */ MPI_Cartdim_get(comm, &cart_dim);
      return cart_dim == 2;
    } else {
      return false;
    }
  }

  void set_blacs_grid() {
    blacs_handle = blacs::sys2blacs_handle(comm);
    blacs_context = blacs_handle;
    char char_grid_major = is_row ? 'R' : 'C';
    blacs::gridinit(blacs_context, char_grid_major, get_nprow(), get_npcol());
  }

private:
  MPI_Comm comm;
  int nprocs, myrank;
  std::array<int,2> size;
  std::array<int,2> my_coordinate;
  bool is_row;
  int blacs_handle, blacs_context;
};

} // namespace rokko

#endif // ROKKO_GRID_HPP
