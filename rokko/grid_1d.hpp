/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2013 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_GRID_1D_HPP
#define ROKKO_GRID_1D_HPP

#include <mpi.h>

namespace rokko {

class grid_1d {
public:
  explicit grid_1d(MPI_Comm comm_in = MPI_COMM_WORLD) : comm(comm_in) {
    initialize();
  }
  MPI_Comm get_comm() const { return comm; }
  int get_nprocs() const { return nprocs; }
  int get_myrank() const { return myrank; }
protected:
  void initialize(void) {
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myrank);
  }
private:
  MPI_Comm comm;
  int nprocs, myrank;
};

} // namespace rokko

#endif // ROKKO_GRID_1D_HPP
