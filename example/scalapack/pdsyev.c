/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
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

  int nb = 1;
  int desc[9];
  int info = cscalapack_descinit(desc, n, n, nb, nb, 0, 0, context, n/nprow);
  if (info) {
    fprintf(stderr, "Error in cscalapack_descinit\n");
    MPI_Abort(MPI_COMM_WORLD, info);
  }
  double **a = alloc_dmatrix(n/nprow, n/npcol);
  double **z = alloc_dmatrix(n/nprow, n/npcol);
  double *w = alloc_dvector(n);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      cscalapack_pdelset(mat_ptr(a), i, j, desc, imin(i, j) + 1);

  info = cscalapack_pdsyev('V', 'U', n, mat_ptr(a), 0, 0, desc, vec_ptr(w),
                           mat_ptr(z), 0, 0, desc);
  if (info) {
    fprintf(stderr, "Error in cscalapack_pdsyev\n");
    MPI_Abort(MPI_COMM_WORLD, info);
  }
  if (myrank == 0) {
    printf("eigenvalues: ");
    fprint_dvector(stdout, n, vec_ptr(w));
  }

  free_dmatrix(a);
  free_dmatrix(z);
  free_dvector(w);
  Cblacs_gridexit(context);
  MPI_Finalize();
}
