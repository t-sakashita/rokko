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

#include <array>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <rokko/mpi_communicator.hpp>
#include <rokko/blacs_grid.hpp>

namespace rokko {

extern struct grid_row_major_t {} grid_row_major;

extern struct grid_col_major_t {} grid_col_major;

class grid : public mpi_comm, public blacs_grid {
public:
  template <typename GRID_MAJOR>
  grid(MPI_Comm comm_in, GRID_MAJOR const& grid_major = grid_row_major)
    : mpi_comm(comm_in) {
    if (comm == MPI_COMM_NULL) {
      set_size({0, 0});  // to avoid zero dividing and blacs does not take (0,0) size grid
      set_my_coordinate({-1,-1});
      return;
    }

    if (is_MPI_2dim_Cart(comm))
      set_sizes_cart(grid_major);
    else
      set_sizes_default(grid_major);

    set_blacs_grid(get_comm(), is_row_major(), get_size());
  }

  template <typename GRID_MAJOR>
  grid(MPI_Comm comm_in, std::array<int,2> const& size_in, GRID_MAJOR const& grid_major = grid_row_major) : mpi_comm(comm_in) {
    if ((std::get<0>(size_in) * std::get<1>(size_in)) != nprocs) {
      throw std::invalid_argument("grid::grid() : (rows * cols) != nprocs");
    }

    set_size(size_in);

    set_major<GRID_MAJOR>();

    set_my_coordinate();

    set_blacs_grid(get_comm(), is_row_major(), get_size());
  }

  explicit grid(MPI_Comm comm_in = MPI_COMM_WORLD) : grid(comm_in, grid_row_major) {}

  static int find_square_root_like_divisor(int n) {
    auto i = static_cast<int>(std::sqrt(static_cast<double>(n)));
    for (; i > 1; --i) {
      if ( (n % i) == 0 ) break;
    }
    return i;
  }

  int get_nprow() const { return std::get<0>(size); }
  int get_npcol() const { return std::get<1>(size); }
  std::array<int,2> get_size() const { return size; }
  int get_myrow() const { return std::get<0>(my_coordinate); }
  int get_mycol() const { return std::get<1>(my_coordinate); }
  std::array<int,2> get_my_coordinate() const { return my_coordinate; }

  void set_size(std::array<int,2> const& size_in) { size = size_in; }
  void set_my_coordinate(std::array<int,2> const& my_coordinate_in) { my_coordinate = my_coordinate_in; }
  void set_my_coordinate() { set_my_coordinate(calculate_coordinate(myrank)); }

  template <typename GRID_MAJOR>
  void set_major() { is_row = std::is_same_v<GRID_MAJOR, grid_row_major_t>; }

  bool is_row_major() const { return is_row; }
  bool is_col_major() const { return !is_row; }

  int calculate_grid_row(int proc_rank) const { 
    return is_row ? proc_rank / size[1] : proc_rank % size[0];
  }
  int calculate_grid_col(int proc_rank) const { 
    return is_row ? proc_rank % size[1] : proc_rank / size[0];
  }

  std::array<int,2> calculate_coordinate(int proc_rank) const {
    return {calculate_grid_row(proc_rank), calculate_grid_col(proc_rank)};
  }

  int calculate_rank_form_coords(int proc_row, int proc_col) const {
    return is_row ? (proc_row * size[1] + proc_col)
      : (proc_col * size[0] + proc_row);
  }

protected:
  template <typename GRID_MAJOR>
  void set_sizes_default(GRID_MAJOR) {
    const auto nprow = find_square_root_like_divisor(nprocs);
    set_size({nprow, nprocs/nprow});
    set_major<GRID_MAJOR>();
    set_my_coordinate();
  }

  static std::string mpi_error_string(int ierr) {
    if (ierr == MPI_ERR_TOPOLOGY)
      return "MPI_ERR_TOPOLOGY";
    else if (ierr == MPI_ERR_COMM)
      return "MPI_ERR_COMM";
    else if (ierr == MPI_ERR_ARG)
      return "MPI_ERR_ARG";
    else
      return "MPI_SUCCESS";
  }

  void set_sizes_cart(grid_row_major_t) {
    constexpr int cart_dim = 2;
    std::array<int,2> dims, periods, coords;
    const auto ierr = MPI_Cart_get(comm, cart_dim, dims.data(), periods.data(), coords.data());
    if (ierr != MPI_SUCCESS)
      throw std::invalid_argument("set_sizes_cart : MPI_Cart_get returns " + mpi_error_string(ierr));

    set_size(dims);
    set_my_coordinate(coords);

    set_major<grid_row_major_t>();
  }

  static void set_sizes_cart(grid_col_major_t) {
    throw std::invalid_argument("calculate_sizes_cart : MPI Cartesian doesn't support grid_col_major.  Use it with grid_row_major.");
  }

private:
  std::array<int,2> size, my_coordinate;
  bool is_row;
};

} // namespace rokko
