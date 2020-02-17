/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_VARIOUS_MPI_COMM_HPP
#define ROKKO_VARIOUS_MPI_COMM_HPP

#include <mpi.h>
#include <array>
#include <rokko/grid.hpp>

MPI_Comm create_even_comm(MPI_Comm comm_in = MPI_COMM_WORLD) {
  MPI_Group group_world, even_group;
  MPI_Comm even_comm;
  int nprocs;

  MPI_Comm_size(comm_in, &nprocs);
  MPI_Comm_group(comm_in, &group_world);

  const int Neven = (nprocs + 1) / 2;
  std::vector<int> even_members(Neven);
  for (int i=0; i<Neven; ++i)
    even_members[i] = 2 * i;

  MPI_Group_incl(group_world, Neven, even_members.data(), &even_group);
  MPI_Comm_create(comm_in, even_group, &even_comm);
  MPI_Group_free(&group_world);
  MPI_Group_free(&even_group);

  return even_comm;
}

MPI_Comm create_odd_comm(MPI_Comm comm_in = MPI_COMM_WORLD) {
  MPI_Group group_world, odd_group;
  MPI_Comm odd_comm;
  int nprocs;

  MPI_Comm_size(comm_in, &nprocs);
  MPI_Comm_group(comm_in, &group_world);

  const int Neven = (nprocs + 1) / 2;
  const int Nodd = nprocs - Neven;
  std::vector<int> odd_members(Nodd);
  for (int i=0; i<Nodd; ++i)
    odd_members[i] = 2 * i + 1;

  MPI_Group_incl(group_world, Nodd, odd_members.data(), &odd_group);
  MPI_Comm_create(comm_in, odd_group, &odd_comm);
  MPI_Group_free(&group_world);
  MPI_Group_free(&odd_group);

  return odd_comm;
}

MPI_Comm create_even_odd_comm(MPI_Comm comm_in = MPI_COMM_WORLD) {
  int rank;
  MPI_Comm_rank(comm_in, &rank);
  MPI_Comm comm = ((rank % 2) == 0) ? create_even_comm(comm_in) : create_odd_comm(comm_in);
  return comm;
}

MPI_Comm create_even_odd_comm_by_split(MPI_Comm comm_in = MPI_COMM_WORLD) {
  MPI_Comm comm;
  int rank;
  MPI_Comm_rank(comm_in, &rank);
  int color = rank % 2;
  MPI_Comm_split(comm_in, color, rank, &comm);
  return comm;
}

std::array<int,2> get_square_like_dims(int nprocs) {
  std::array<int,2> dims;
  dims[0] = rokko::grid::find_square_root_like_divisor(nprocs);
  dims[1] = nprocs / dims[0];

  return dims;
}

MPI_Comm create_cart_comm(MPI_Comm comm_in = MPI_COMM_WORLD) {
  int nprocs;
  MPI_Comm_size(comm_in, &nprocs);

  MPI_Comm comm;
  std::array<int,2> dims = get_square_like_dims(nprocs);
  std::array<int,2> periods{0,0};
  int reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims.data(), periods.data(), reorder, &comm);

  return comm;
}

#endif // ROKKO_VARIOUS_MPI_COMM_HPP
