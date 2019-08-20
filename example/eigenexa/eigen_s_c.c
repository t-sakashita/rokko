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
#include <stdio.h>
#include <rokko/ceigenexa.h>
#include <rokko/cmatrix.h>

int imin(int x, int y) { return (x < y) ? x : y; }

int main(int argc, char** argv) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  ceigenexa_init();

  int n = 8;
  if (argc > 1) n = atoi(argv[1]);

  int nprocs, npcol, nprow;
  ceigenexa_get_procs(&nprocs, &npcol, &nprow);
  if (myrank == 0) {
    printf("n = %d\n", n);
    printf("nprocs = %d\n", nprocs);
    printf("nprow = %d\n", nprow);
    printf("npcol = %d\n", npcol);
  }
  int inod, x_inod, y_inod;
  ceigenexa_get_id(&inod, &x_inod, &y_inod);

  int nm, ny;
  ceigenexa_get_matdims(nprow, npcol, n, &nm, &ny);
  double **a = alloc_dmatrix(nm, ny);
  double **z = alloc_dmatrix(nm, ny);
  double *w = alloc_dvector(n);

  int jloop_sta = ceigenexa_loop_start(0, nprow, y_inod);
  int jloop_end = ceigenexa_loop_end(n, nprow, y_inod);
  int iloop_sta = ceigenexa_loop_start(0, npcol, x_inod);
  int iloop_end = ceigenexa_loop_end(n, npcol, x_inod);
  for (int jl = jloop_sta; jl < jloop_end; ++jl) {
    int j = ceigenexa_translate_l2g(jl, nprow, y_inod);
    for (int il = iloop_sta; il < iloop_end; ++il) {
      int i = ceigenexa_translate_l2g(il, npcol, x_inod);
      mat_elem(a, il, jl) = imin(i, j) + 1;
    }
  }

  ceigenexa_eigen_s(n, n, mat_ptr(a), nm, vec_ptr(w), mat_ptr(z), nm, 48, 128, 'A');
  if (myrank == 0) {
    printf("eigenvalues: ");
    fprint_dvector(stdout, n, vec_ptr(w));
  }

  free_dmatrix(a);
  free_dmatrix(z);
  free_dvector(w);
  ceigenexa_free();
  MPI_Finalize();
}
