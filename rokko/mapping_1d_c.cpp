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

#include <rokko/mapping_1d.h>
#include <rokko/mpi_communicator.hpp>
#include <rokko/mapping_1d.hpp>

void rokko_mapping_1d_construct(struct rokko_mapping_1d* map, int dim, MPI_Comm comm) {
  map->ptr = new rokko::mapping_1d(dim, rokko::mpi_comm{comm});
}

void rokko_mapping_1d_construct_f(struct rokko_mapping_1d* map, int dim, int comm_f) {
  MPI_Comm comm = MPI_Comm_f2c(comm_f);
  map->ptr = new rokko::mapping_1d(dim, rokko::mpi_comm{comm});
}

void rokko_mapping_1d_destruct(struct rokko_mapping_1d* map) {
  delete static_cast<rokko::mapping_1d*>(map->ptr);
  map->ptr = nullptr;
}

int rokko_mapping_1d_get_dim(struct rokko_mapping_1d map) {
  return static_cast<rokko::mapping_1d*>(map.ptr)->get_dim();
}

int rokko_mapping_1d_num_local_rows(struct rokko_mapping_1d map) {
  return static_cast<rokko::mapping_1d*>(map.ptr)->num_local_rows();
}

int rokko_mapping_1d_start_row(struct rokko_mapping_1d map) {
  return static_cast<rokko::mapping_1d*>(map.ptr)->start_row();
}

int rokko_mapping_1d_end_row(struct rokko_mapping_1d map) {
  return static_cast<rokko::mapping_1d*>(map.ptr)->end_row();
}

MPI_Comm rokko_mapping_1d_get_comm(struct rokko_mapping_1d map) {
  return static_cast<rokko::mapping_1d*>(map.ptr)->get_mpi_comm().get_comm();
}

MPI_Fint rokko_mapping_1d_get_comm_f(struct rokko_mapping_1d map) {
  return MPI_Comm_c2f(rokko_mapping_1d_get_comm(map));
}
