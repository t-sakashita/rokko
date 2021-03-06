/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rokko/utility/heisenberg_hamiltonian_mpi.h>

#define min(a, b) ((a) < (b) ? (a) : (b))

struct heisenberg_vars {
  MPI_Comm comm;
  int nprocs, myrank;
  int dim;
  int L;
  int lattice_size;
  int *lattice_first, *lattice_second;
  double *buffer;
};

// Assume that n >= 1
int find_power_of_two(int n) {
  int p = 0;
  while (n > 1) {
    n /= 2;
    ++p;
  };
  return p;
}

void heisenberg_multiply(const double *const x, double *const y, void* vars) {
  struct heisenberg_vars* p = (struct heisenberg_vars*)vars;
  multiply(p->comm, p->lattice_size, p->L, p->lattice_first, p->lattice_second, x, y, p->buffer);
}

void heisenberg_initialize(struct rokko_distributed_mfree* mat, int L, int lattice_size, int lattice_first[], int lattice_second[], struct heisenberg_vars* vars) {
  vars->comm = MPI_COMM_WORLD;
  MPI_Comm_size(vars->comm, &vars->nprocs);
  MPI_Comm_rank(vars->comm, &vars->myrank);

  vars->L = L;
  vars->lattice_size = lattice_size;
  vars->lattice_first = lattice_first;
  vars->lattice_second = lattice_second;

  int dim = 1 << L;
  rokko_distributed_mfree_construct(mat, heisenberg_multiply, vars, dim, MPI_COMM_WORLD);

  int num_local_rows = rokko_distributed_mfree_num_local_rows(*mat);
  vars->buffer = (double*) malloc(sizeof(double) * num_local_rows);
}


int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char* library_routine, *library, *routine;
  if (argc >= 2) {
    library_routine = rokko_parallel_sparse_ev_default_solver();
    if (argc >= 2) library_routine = argv[1];
  } else {
    library_routine = rokko_parallel_sparse_ev_default_solver();
  }
  rokko_split_solver_name(library_routine, &library, &routine);
  int L = (argc == 3) ? L = atoi(argv[2]) : 10;
  int lattice_size = L;
  int* lattice_first = (int*) malloc(sizeof(int) * lattice_size);
  int* lattice_second = (int*) malloc(sizeof(int) * lattice_size);
  int i;
  for (i = 0; i < L; ++i) {
    lattice_first[i] = i;
    lattice_second[i] = (i+1) % L;
  }

  struct rokko_parallel_sparse_ev solver;
  rokko_parallel_sparse_ev_construct(&solver, library, argc, argv);
  struct rokko_distributed_mfree mat;
  struct heisenberg_vars vars;
  heisenberg_initialize(&mat, L, lattice_size, lattice_first, lattice_second, &vars);
  int dim = rokko_distributed_mfree_dim(mat);
  if (rank == 0) {
    printf("Eigenvalue decomposition of antiferromagnetic Heisenberg chain\n");
    printf("solver = %s\n", library);
    printf("dimension = %d\n", dim);
  }

  struct rokko_parameters params;
  rokko_parameters_construct(&params);
  // set some parameters
  if (routine[0] != '\0')  rokko_parameters_set_string(params, "routine", routine);
  rokko_parameters_set_int(params, "block_size", 5);
  rokko_parameters_set_int(params, "max_iters", 500);
  rokko_parameters_set_double(params, "conv_tol", 1.0e-8);
  rokko_parameters_set_int(params, "num_eigenvalues", 10);
  rokko_parallel_sparse_ev_diagonalize_distributed_mfree(solver, mat, params);

  int num_conv = rokko_parallel_sparse_ev_num_conv(solver);
  if (num_conv == 0) MPI_Abort(MPI_COMM_WORLD, -1);
  int num_local_rows = rokko_distributed_mfree_num_local_rows(mat);
  double eig_vec[num_local_rows];
  rokko_parallel_sparse_ev_eigenvector(solver, 0, eig_vec);
  if (rank == 0) {
    printf("number of converged eigenpairs = %d\n", num_conv);
    printf("smallest eigenvalues: ");
    int i, j;
    for (i = 0; i < num_conv; ++i) printf("%g", rokko_parallel_sparse_ev_eigenvalue(solver, i));
    printf("\n");
    printf("smallest eigenvector: ");
    for (j = 0; j < num_local_rows; ++j)
      printf("%.5e ", eig_vec[j]);
    printf("\n");
  }
  rokko_distributed_mfree_destruct(&mat);
  rokko_parallel_sparse_ev_destruct(&solver);

  MPI_Finalize();
  return 0;
}
