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

#include <rokko/parallel_sparse_solver.hpp>
#include <rokko/rokko_sparse.h>
#include <rokko/parameters.h>

void rokko_parallel_sparse_solver_construct(struct rokko_parallel_sparse_solver* solver, const char* solver_name, int argc, char** argv) {
  solver->ptr = new rokko::parallel_sparse_solver(std::string(solver_name));
  static_cast<rokko::parallel_sparse_solver*>(solver->ptr)->initialize(argc, argv);
}

void rokko_parallel_sparse_solver_construct_f(rokko_parallel_sparse_solver* solver, const char* solver_name) {
  int argc = 0;
  char** argv;
  rokko_parallel_sparse_solver_construct(solver, solver_name, argc, argv);
}

void rokko_parallel_sparse_solver_destruct(rokko_parallel_sparse_solver* solver) {
  rokko::parallel_sparse_solver* ptr = static_cast<rokko::parallel_sparse_solver*>(solver->ptr);
  ptr->finalize();
  delete ptr;
}

struct rokko_parameters rokko_parallel_sparse_solver_diagonalize_distributed_crs_matrix(struct rokko_parallel_sparse_solver* solver,
								     struct rokko_distributed_crs_matrix* mat,
								     int num_evals, int block_size, int max_iters, double tol) {
  struct rokko_parameters params_out;
  rokko_parameters_construct(&params_out);
  *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::parallel_sparse_solver*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_crs_matrix*>(mat->ptr),
									num_evals, block_size, max_iters, tol);
  return params_out;
}

struct rokko_parameters rokko_parallel_sparse_solver_diagonalize_distributed_mfree(struct rokko_parallel_sparse_solver* solver,
								struct rokko_distributed_mfree* mat,
								int num_evals, int block_size, int max_iters, double tol) {
  struct rokko_parameters params_out;
  rokko_parameters_construct(&params_out);
  *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::parallel_sparse_solver*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_mfree*>(mat->ptr),
									num_evals, block_size, max_iters, tol);
  return params_out;
}

double rokko_parallel_sparse_solver_eigenvalue(struct rokko_parallel_sparse_solver* solver, int i) {
  return static_cast<rokko::parallel_sparse_solver*>(solver->ptr)->eigenvalue(i);
}

void rokko_parallel_sparse_solver_eigenvector(struct rokko_parallel_sparse_solver* solver, int i, double* vec) {
  return static_cast<rokko::parallel_sparse_solver*>(solver->ptr)->eigenvector(i, vec);
}

int rokko_parallel_sparse_solver_num_conv(struct rokko_parallel_sparse_solver* solver) {
  return static_cast<rokko::parallel_sparse_solver*>(solver->ptr)->num_conv();
}

int rokko_parallel_sparse_solver_num_solvers() {
  return rokko::parallel_sparse_solver::solvers().size();
}

char** rokko_parallel_sparse_solver_solvers() {
  std::vector<std::string> solvers = rokko::parallel_sparse_solver::solvers();
  char **solvers_c = (char**)malloc((size_t)(solvers.size() * sizeof(char*)));
  for (int i = 0; i < solvers.size(); ++i) {
    solvers_c[i] = (char*)malloc((size_t)((solvers[i].size() + 1) * sizeof(char)));
    strcpy(solvers_c[i], solvers[i].c_str());
  }
  return solvers_c;
}

char* rokko_parallel_sparse_solver_default_solver() {
  std::string solver = rokko::parallel_sparse_solver::default_solver();
  char *solver_c = (char*)malloc((size_t)((solver.size() + 1) * sizeof(char)));
  strcpy(solver_c, solver.c_str());
  return solver_c;
}
