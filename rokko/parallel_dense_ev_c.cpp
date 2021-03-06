/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/parallel_dense_ev.h>
#include <rokko/parallel_dense_ev.hpp>
#include <rokko/copy_string.hpp>

void rokko_parallel_dense_ev_construct(struct rokko_parallel_dense_ev* solver, const char* solver_name, int argc, char** argv) {
  solver->ptr = new rokko::parallel_dense_ev(std::string(solver_name));
  static_cast<rokko::parallel_dense_ev*>(solver->ptr)->initialize(argc, argv);
}

void rokko_parallel_dense_ev_construct_f(struct rokko_parallel_dense_ev* solver, const char* solver_name) {
  int argc = 0;
  char** argv = 0;
  rokko_parallel_dense_ev_construct(solver, solver_name, argc, argv);
}

void rokko_parallel_dense_ev_destruct(struct rokko_parallel_dense_ev* solver) {
  rokko::parallel_dense_ev* ptr = static_cast<rokko::parallel_dense_ev*>(solver->ptr);
  ptr->finalize();
  delete ptr;
  solver->ptr = nullptr;
}

struct rokko_parameters rokko_parallel_dense_ev_diagonalize(struct rokko_parallel_dense_ev solver,
  struct rokko_distributed_matrix mat, struct rokko_eigen_vector eigvals,
  struct rokko_distributed_matrix eigvecs, struct rokko_parameters params) {
  struct rokko_parameters params_out;
  rokko_parameters_construct(&params_out);

  if (mat.major == rokko_matrix_col_major)
    *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::parallel_dense_ev*>(solver.ptr)->diagonalize(
      *static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(mat.ptr),
      *static_cast<Eigen::VectorXd*>(eigvals.ptr),
      *static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(eigvecs.ptr),
      *static_cast<rokko::parameters*>(params.ptr));
  else
    *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::parallel_dense_ev*>(solver.ptr)->diagonalize(
      *static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(mat.ptr),
      *static_cast<Eigen::VectorXd*>(eigvals.ptr),
      *static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(eigvecs.ptr),
      *static_cast<rokko::parameters*>(params.ptr));
  return params_out;
}

struct rokko_parameters rokko_parallel_dense_ev_diagonalize_distributed_matrix(struct rokko_parallel_dense_ev solver,
  struct rokko_distributed_matrix mat, struct rokko_eigen_vector eigvals,
  struct rokko_distributed_matrix eigvecs) {
  struct rokko_parameters params_out;
  rokko_parameters_construct(&params_out);

  if (mat.major == rokko_matrix_col_major)
    *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::parallel_dense_ev*>(solver.ptr)->diagonalize(
      *static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(mat.ptr),
      *static_cast<Eigen::VectorXd*>(eigvals.ptr),
      *static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(eigvecs.ptr));
  else
    *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::parallel_dense_ev*>(solver.ptr)->diagonalize(
      *static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(mat.ptr),
      *static_cast<Eigen::VectorXd*>(eigvals.ptr),
      *static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(eigvecs.ptr));
  return params_out;
}

// For Fortran binding
// eigenvalues/eigenvectors, parameters
void rokko_parallel_dense_ev_diagonalize_f(struct rokko_parallel_dense_ev* solver,
					   struct rokko_distributed_matrix* mat, struct rokko_eigen_vector* eigvals,
					   struct rokko_distributed_matrix* eigvecs,
					   struct rokko_parameters* params, struct rokko_parameters* params_out) {
  rokko_parameters_construct(params_out);
  if (mat->major == rokko_matrix_col_major)
    *static_cast<rokko::parameters*>(params_out->ptr) = static_cast<rokko::parallel_dense_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(mat->ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								       *static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(eigvecs->ptr),
								       *static_cast<rokko::parameters*>(params->ptr));
  else
    *static_cast<rokko::parameters*>(params_out->ptr) = static_cast<rokko::parallel_dense_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(mat->ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								       *static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(eigvecs->ptr),
								       *static_cast<rokko::parameters*>(params->ptr));
}

void rokko_parallel_dense_ev_diagonalize_no_params_out_f(struct rokko_parallel_dense_ev* solver,
						       struct rokko_distributed_matrix* mat, struct rokko_eigen_vector* eigvals,
						       struct rokko_distributed_matrix* eigvecs,
						       struct rokko_parameters* params) {
  if (mat->major == rokko_matrix_col_major)
    static_cast<rokko::parallel_dense_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								   *static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(eigvecs->ptr),
								   *static_cast<rokko::parameters*>(params->ptr));
  else
    static_cast<rokko::parallel_dense_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								   *static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(eigvecs->ptr),
								   *static_cast<rokko::parameters*>(params->ptr));
}


void rokko_parallel_dense_ev_diagonalize_no_params_inout_f(struct rokko_parallel_dense_ev* solver,
							 struct rokko_distributed_matrix* mat, struct rokko_eigen_vector* eigvals,
							 struct rokko_distributed_matrix* eigvecs) {
  if (mat->major == rokko_matrix_col_major)
    static_cast<rokko::parallel_dense_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								   *static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(eigvecs->ptr));
  else
    static_cast<rokko::parallel_dense_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								   *static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(eigvecs->ptr));
}

void rokko_parallel_dense_ev_diagonalize_eigvals_f(struct rokko_parallel_dense_ev* solver,
						 struct rokko_distributed_matrix* mat, struct rokko_eigen_vector* eigvals,
						 struct rokko_parameters* params, struct rokko_parameters* params_out) {
  rokko_parameters_construct(params_out);
  if (mat->major == rokko_matrix_col_major)
    *static_cast<rokko::parameters*>(params_out->ptr) = static_cast<rokko::parallel_dense_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(mat->ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								       *static_cast<rokko::parameters*>(params->ptr));
  else
    *static_cast<rokko::parameters*>(params_out->ptr) = static_cast<rokko::parallel_dense_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(mat->ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								       *static_cast<rokko::parameters*>(params->ptr));
}

void rokko_parallel_dense_ev_diagonalize_eigvals_no_params_out_f(struct rokko_parallel_dense_ev* solver,
							       struct rokko_distributed_matrix* mat, struct rokko_eigen_vector* eigvals,
							       struct rokko_parameters* params) {
  if (mat->major == rokko_matrix_col_major)
    static_cast<rokko::parallel_dense_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								   *static_cast<rokko::parameters*>(params->ptr));
  else
    static_cast<rokko::parallel_dense_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								   *static_cast<rokko::parameters*>(params->ptr));
}

void rokko_parallel_dense_ev_diagonalize_eigvals_no_params_inout_f(struct rokko_parallel_dense_ev* solver,
								 struct rokko_distributed_matrix* mat, struct rokko_eigen_vector* eigvals) {
  if (mat->major == rokko_matrix_col_major)
    static_cast<rokko::parallel_dense_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr));
  else
    static_cast<rokko::parallel_dense_ev*>(solver->ptr)->diagonalize(*static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr));
}

struct rokko_mapping_bc rokko_parallel_dense_ev_default_mapping(struct rokko_parallel_dense_ev solver, int dim, struct rokko_grid g) {
  struct rokko_mapping_bc map;
  map.major = rokko_matrix_col_major;
  map.ptr = new rokko::mapping_bc<rokko::matrix_col_major>(static_cast<rokko::parallel_dense_ev*>(solver.ptr)->default_mapping(dim, *static_cast<rokko::grid*>(g.ptr)));
  return map;
}

void rokko_parallel_dense_ev_default_mapping_f(struct rokko_parallel_dense_ev solver, int dim, struct rokko_grid g, struct rokko_mapping_bc* map) {
  map->major = rokko_matrix_col_major;
  map->ptr = new rokko::mapping_bc<rokko::matrix_col_major>(static_cast<rokko::parallel_dense_ev*>(solver.ptr)->default_mapping(dim, *static_cast<rokko::grid*>(g.ptr)));
}

int rokko_parallel_dense_ev_num_solvers() {
  return rokko::parallel_dense_ev::solvers().size();
}

char** rokko_parallel_dense_ev_solvers() {
  std::vector<std::string> solvers = rokko::parallel_dense_ev::solvers();
  char **solvers_c = (char**)malloc((size_t)(solvers.size() * sizeof(char*)));
  for (std::size_t i = 0; i < solvers.size(); ++i) {
    solvers_c[i] = copy_string(solvers[i]);
  }
  return solvers_c;
}

const char* rokko_parallel_dense_ev_default_solver() {
  return rokko::parallel_dense_ev::default_solver().c_str();
}
