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

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <rokko/Cblacs.h>
#include <rokko/cscalapack.h>
#include <rokko/cmatrix.h>

int imin(int x, int y) { return (x < y) ? x : y; }

int main(int argc, char** argv) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int nprocs, myrank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int n = 8;
  if (argc > 1) n = atoi(argv[1]);

  int nprow = (int)(sqrt(1. * nprocs + 0.5));
  int npcol = (int)(1. * nprocs / nprow);
  if (myrank == 0) {
    printf("n = %d\n", n);
    printf("nprocs = %d\n", nprocs);
    printf("nprow = %d\n", nprow);
    printf("npcol = %d\n", npcol);
    if ((nprocs != nprow * npcol) || (n % nprow != 0) || (n % npcol != 0)) {
      printf("incompatible matrix size and number of processes");
      MPI_Abort(MPI_COMM_WORLD, 127);
    }
  }
  int context;
  Cblacs_get(0, 0, &context);
  char order = 'R';
  Cblacs_gridinit(&context, &order, nprow, npcol);
  int myrow, mycol;
  Cblacs_gridinfo(context, &nprow, &npcol, &myrow, &mycol);

  int nb = 1;
  int desc[9];
  int info = cscalapack_descinit(desc, n, n, nb, nb, 0, 0, context, n/nprow);
  if (info) {
    fprintf(stderr, "Error in cscalapack_descinit\n");
    MPI_Abort(MPI_COMM_WORLD, info);
  }
  double **a = alloc_dmatrix(n/nprow, n/npcol);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      cscalapack_pdelset(mat_ptr(a), i, j, desc, imin(i, j) + 1);

  double work[1000];
  const char cmatnm[] = "minij_matrix";
  cscalapack_pdlaprnt(n, n, mat_ptr(a), 0, 0, desc, myrow, mycol, cmatnm, 6, work);

  free_dmatrix(a);
  Cblacs_gridexit(context);
  MPI_Finalize();
}
