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
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int color = rank % 2;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);
  return comm;
}

#endif // ROKKO_VARIOUS_MPI_COMM_HPP
