/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/ceigenexa.h>
#include <rokko/eigenexa/eigenexa_interface.h>

void ceigenexa_init(void) {
  EIGENEXA_init_wrap0();
}

void ceigenexa_init1(MPI_Comm comm) {
  MPI_Fint fcomm = MPI_Comm_c2f(comm);
  EIGENEXA_init_wrap1(&fcomm);
}

void ceigenexa_init2(MPI_Comm comm, char grid_major) {
  MPI_Fint fcomm = MPI_Comm_c2f(comm);
  EIGENEXA_init_wrap2(&fcomm, &grid_major);
}
