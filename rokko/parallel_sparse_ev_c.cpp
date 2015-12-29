/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>,
*                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/parallel_sparse_ev.hpp>
#include <rokko/rokko_sparse.h>
#include <rokko/parameters.h>

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
  solver->ptr = 0;
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
  *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::parallel_sparse_ev*>(solver.ptr)->diagonalize(*static_cast<rokko::distributed_mfree*>(mat.ptr),
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
  *static_cast<rokko::parameters*>(params_out->ptr) = static_cast<rokko::parallel_sparse_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_mfree*>(mat->ptr),
															*static_cast<rokko::parameters*>(params->ptr));
}

void rokko_parallel_sparse_ev_diagonalize_distributed_mfree_noreturn_f(struct rokko_parallel_sparse_ev* solver,
								       struct rokko_distributed_mfree* mat,
								       struct rokko_parameters* params) {
  static_cast<rokko::parallel_sparse_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_mfree*>(mat->ptr),
								    *static_cast<rokko::parameters*>(params->ptr));
}


double rokko_parallel_sparse_ev_eigenvalue(struct rokko_parallel_sparse_ev solver, int i) {
  return static_cast<rokko::parallel_sparse_ev*>(solver.ptr)->eigenvalue(i);
}

void rokko_parallel_sparse_ev_eigenvector(struct rokko_parallel_sparse_ev solver, int i, double* vec) {
  return static_cast<rokko::parallel_sparse_ev*>(solver.ptr)->eigenvector(i, vec);
}

int rokko_parallel_sparse_ev_num_conv(struct rokko_parallel_sparse_ev solver) {
  return static_cast<rokko::parallel_sparse_ev*>(solver.ptr)->num_conv();
}

int rokko_parallel_sparse_ev_num_solvers() {
  return rokko::parallel_sparse_ev::solvers().size();
}

char** rokko_parallel_sparse_ev_solvers() {
  std::vector<std::string> solvers = rokko::parallel_sparse_ev::solvers();
  char **solvers_c = (char**)malloc((size_t)(solvers.size() * sizeof(char*)));
  for (int i = 0; i < solvers.size(); ++i) {
    solvers_c[i] = (char*)malloc((size_t)((solvers[i].size() + 1) * sizeof(char)));
    strcpy(solvers_c[i], solvers[i].c_str());
  }
  return solvers_c;
}

char* rokko_parallel_sparse_ev_default_solver() {
  std::string solver = rokko::parallel_sparse_ev::default_solver();
  char *solver_c = (char*)malloc((size_t)((solver.size() + 1) * sizeof(char)));
  strcpy(solver_c, solver.c_str());
  return solver_c;
}
