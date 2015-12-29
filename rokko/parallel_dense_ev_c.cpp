/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/solver.hpp>
#include <rokko/rokko_dense.h>
#include <rokko/parameters.h>

void rokko_parallel_dense_ev_construct(rokko_parallel_dense_ev* solver, const char* solver_name, int argc, char** argv) {
  solver->ptr = new rokko::parallel_dense_ev(std::string(solver_name));
  static_cast<rokko::parallel_dense_ev*>(solver->ptr)->initialize(argc, argv);
}

void rokko_parallel_dense_ev_construct_f(rokko_parallel_dense_ev* solver, const char* solver_name) {
  int argc = 0;
  char** argv;
  rokko_parallel_dense_ev_construct(solver, solver_name, argc, argv);
}

void rokko_parallel_dense_ev_destruct(rokko_parallel_dense_ev* solver) {
  rokko::parallel_dense_ev* ptr = static_cast<rokko::parallel_dense_ev*>(solver->ptr);
  ptr->finalize();
  delete ptr;
}

struct rokko_parameters rokko_parallel_dense_ev_diagonalize_distributed_matrix(struct rokko_parallel_dense_ev solver,
  struct rokko_distributed_matrix mat, struct rokko_localized_vector eigvals,
  struct rokko_distributed_matrix eigvecs) {
  struct rokko_parameters params_out;
  rokko_parameters_construct(&params_out);

  if (mat.major == rokko_matrix_col_major)
    *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::parallel_dense_ev*>(solver.ptr)->diagonalize(
      *static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(mat.ptr),
      *static_cast<rokko::localized_vector<double>*>(eigvals.ptr),
      *static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(eigvecs.ptr));
  else
    *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::parallel_dense_ev*>(solver.ptr)->diagonalize(
      *static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(mat.ptr),
      *static_cast<rokko::localized_vector<double>*>(eigvals.ptr),
      *static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(eigvecs.ptr));
  return params_out;
}

int rokko_parallel_dense_ev_num_solvers() {
  return rokko::parallel_dense_ev::solvers().size();
}

char** rokko_parallel_dense_ev_solvers() {
  std::vector<std::string> solvers = rokko::parallel_dense_ev::solvers();
  char **solvers_c = (char**)malloc((size_t)(solvers.size() * sizeof(char*)));
  for (int i = 0; i < solvers.size(); ++i) {
    solvers_c[i] = (char*)malloc((size_t)((solvers[i].size() + 1) * sizeof(char)));
    strcpy(solvers_c[i], solvers[i].c_str());
  }
  return solvers_c;
}

char* rokko_parallel_dense_ev_default_solver() {
  std::string solver = rokko::parallel_dense_ev::default_solver();
  char *solver_c = (char*)malloc((size_t)((solver.size() + 1) * sizeof(char)));
  strcpy(solver_c, solver.c_str());
  return solver_c;
}
