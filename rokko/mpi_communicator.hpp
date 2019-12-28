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

#ifndef ROKKO_MPI_COMMNUNICATOR_HPP
#define ROKKO_MPI_COMMNUNICATOR_HPP

#include <mpi.h>

namespace rokko {

class mpi_comm {
public:
  explicit mpi_comm(MPI_Comm comm_in = MPI_COMM_WORLD) : comm(comm_in) {
    if (comm != MPI_COMM_NULL) {
      MPI_Comm_size(comm, &nprocs);
      MPI_Comm_rank(comm, &myrank);
    } else {
      nprocs = 0;
      myrank = -1;
    }
  }

  MPI_Comm get_comm() const { return comm; }
  int get_nprocs() const { return nprocs; }
  int get_myrank() const { return myrank; }

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

protected:
  MPI_Comm comm;
  int nprocs, myrank;
};

} // namespace rokko

#endif // ROKKO_MPI_COMMNUNICATOR_HPP
