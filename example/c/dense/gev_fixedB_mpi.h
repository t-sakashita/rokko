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

#ifndef ROKKO_GEV_FIXEDB_MPI_H
#define ROKKO_GEV_FIXEDB_MPI_H

#include <mpi.h>
#include <rokko/rokko.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void diagonalize_fixedB_c(struct rokko_parallel_dense_ev solver_in, struct rokko_distributed_matrix A_in, struct rokko_distributed_matrix B,
			  struct rokko_eigen_vector eigval_in, struct rokko_distributed_matrix eigvec_in, double tol);

void set_A_B_c(struct rokko_eigen_matrix locA_in, struct rokko_eigen_matrix locB_in);

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_GEV_FIXEDB_MPI_H */


