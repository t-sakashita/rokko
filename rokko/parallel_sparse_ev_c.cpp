/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/parallel_sparse_ev.h>
#include <rokko/parallel_sparse_ev.hpp>
#include <rokko/copy_string.hpp>

void rokko_parallel_sparse_ev_construct(struct rokko_parallel_sparse_ev* solver, const char* solver_name, int argc, char** argv) {
  solver->ptr = new rokko::parallel_sparse_ev(std::string(solver_name));
  static_cast<rokko::parallel_sparse_ev*>(solver->ptr)->initialize(argc, argv);
}

void rokko_parallel_sparse_ev_construct_f(rokko_parallel_sparse_ev* solver, const char* solver_name) {
  int argc = 0;
  char** argv;
  rokko_parallel_sparse_ev_construct(solver, solver_name, argc, argv);
}

void rokko_parallel_sparse_ev_destruct(rokko_parallel_sparse_ev* solver) {
  rokko::parallel_sparse_ev* ptr = static_cast<rokko::parallel_sparse_ev*>(solver->ptr);
  ptr->finalize();
  delete ptr;
  solver->ptr = nullptr;
}


struct rokko_parameters rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix(struct rokko_parallel_sparse_ev solver,
										    struct rokko_distributed_crs_matrix mat,
										    struct rokko_parameters params) {
  struct rokko_parameters params_out;
  rokko_parameters_construct(&params_out);
  *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::parallel_sparse_ev*>(solver.ptr)->diagonalize(*static_cast<rokko::distributed_crs_matrix*>(mat.ptr),
														      *static_cast<rokko::parameters*>(params.ptr));
  return params_out;
}

struct rokko_parameters rokko_parallel_sparse_ev_diagonalize_distributed_mfree(struct rokko_parallel_sparse_ev solver,
									       struct rokko_distributed_mfree mat,
									       struct rokko_parameters params) {
  struct rokko_parameters params_out;
  rokko_parameters_construct(&params_out);
  *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::parallel_sparse_ev*>(solver.ptr)->diagonalize(*static_cast<rokko::distributed_mfree_holder*>(mat.ptr),
														       *static_cast<rokko::parameters*>(params.ptr));
  return params_out;
}

// Fortran binding
void rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix_f(struct rokko_parallel_sparse_ev* solver,
								   struct rokko_distributed_crs_matrix* mat,
								   struct rokko_parameters* params, struct rokko_parameters* params_out) {
  rokko_parameters_construct(params_out);
  *static_cast<rokko::parameters*>(params_out->ptr) = static_cast<rokko::parallel_sparse_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_crs_matrix*>(mat->ptr),
															*static_cast<rokko::parameters*>(params->ptr));
}

void rokko_parallel_sparse_ev_diagonalize_distributed_crs_matrix_noreturn_f(struct rokko_parallel_sparse_ev* solver,
									    struct rokko_distributed_crs_matrix* mat,
									    struct rokko_parameters* params) {
  static_cast<rokko::parallel_sparse_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_crs_matrix*>(mat->ptr),
								    *static_cast<rokko::parameters*>(params->ptr));
}

void rokko_parallel_sparse_ev_diagonalize_distributed_mfree_f(struct rokko_parallel_sparse_ev* solver,
							      struct rokko_distributed_mfree* mat,
							      struct rokko_parameters* params, struct rokko_parameters* params_out) {
  rokko_parameters_construct(params_out);
  *static_cast<rokko::parameters*>(params_out->ptr) = static_cast<rokko::parallel_sparse_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_mfree_holder*>(mat->ptr),
															*static_cast<rokko::parameters*>(params->ptr));
}

void rokko_parallel_sparse_ev_diagonalize_distributed_mfree_noreturn_f(struct rokko_parallel_sparse_ev* solver,
								       struct rokko_distributed_mfree* mat,
								       struct rokko_parameters* params) {
  static_cast<rokko::parallel_sparse_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_mfree_holder*>(mat->ptr),
								    *static_cast<rokko::parameters*>(params->ptr));
}


double rokko_parallel_sparse_ev_eigenvalue(struct rokko_parallel_sparse_ev solver, int i) {
  return static_cast<rokko::parallel_sparse_ev*>(solver.ptr)->eigenvalue(i);
}

void rokko_parallel_sparse_ev_eigenvector(struct rokko_parallel_sparse_ev solver, int i, double* vec) {
  return static_cast<rokko::parallel_sparse_ev*>(solver.ptr)->eigenvector(i, vec);
}

int rokko_parallel_sparse_ev_num_conv(struct rokko_parallel_sparse_ev solver) {
  return static_cast<rokko::parallel_sparse_ev*>(solver.ptr)->get_num_conv();
}

struct rokko_mapping_1d rokko_parallel_sparse_ev_default_mapping(struct rokko_parallel_sparse_ev solver, int dim, MPI_Comm comm) {
  struct rokko_mapping_1d map;
  map.ptr = new rokko::mapping_1d(static_cast<rokko::parallel_sparse_ev*>(solver.ptr)->default_mapping(dim, rokko::mpi_comm{comm}));
  return map;
}

void rokko_parallel_sparse_ev_default_mapping_f(struct rokko_parallel_sparse_ev solver, int dim, int comm_f, struct rokko_mapping_1d* map) {
  MPI_Comm comm = MPI_Comm_f2c(comm_f);
  map->ptr = new rokko::mapping_1d(static_cast<rokko::parallel_sparse_ev*>(solver.ptr)->default_mapping(dim, rokko::mpi_comm{comm}));
}

int rokko_parallel_sparse_ev_num_solvers() {
  return rokko::parallel_sparse_ev::solvers().size();
}

char** rokko_parallel_sparse_ev_solvers() {
  std::vector<std::string> solvers = rokko::parallel_sparse_ev::solvers();
  char **solvers_c = (char**)malloc((size_t)(solvers.size() * sizeof(char*)));
  for (int i = 0; i < solvers.size(); ++i) {
    solvers_c[i] = copy_string(solvers[i]);
  }
  return solvers_c;
}

const char* rokko_parallel_sparse_ev_default_solver() {
  return rokko::parallel_sparse_ev::default_solver().c_str();
}
