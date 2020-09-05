/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
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


void laplacian_generate(int n, double* d, double* e) {
  int i;

  d[0] = 1;
  for(i = 1; i < n; ++i)
    d[i] = 2;
  for(i = 0; i < (n-1); ++i)
    e[i] = -1;
}

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

  double *d = alloc_dvector(n);
  double *e = alloc_dvector(n);
  laplacian_generate(n, d, e);

  double *w = alloc_dvector(n);

  char range = 'A', order_stored = 'E';
  double vl = 0, vu = 0.;
  int il = 0, iu = 0;
  double abstol = 0.;
  int m, nsplit;
  int *iblock = alloc_ivector(n);
  int *isplit = alloc_ivector(n);

  int info = cscalapack_pdstebz(context, range, order_stored, n,
                                vl, vu, il, iu,
                                abstol, d, e, &m, &nsplit,
                                w, iblock, isplit);

  if (info) {
    fprintf(stderr, "Error in cscalapack_pdstebz\n");
    MPI_Abort(MPI_COMM_WORLD, info);
  }
  if (myrank == 0) {
    printf("eigenvalues: ");
    fprint_dvector(stdout, n, vec_ptr(w));
  }

  free_ivector(iblock);
  free_ivector(isplit);
  free_dvector(w);
  free_dvector(d);
  free_dvector(e);
  Cblacs_gridexit(context);
  MPI_Finalize();
}
