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

#include <rokko/serial_dense_ev.h>
#include <rokko/serial_dense_ev.hpp>
#include <rokko/copy_string.hpp>

void rokko_serial_dense_ev_construct(rokko_serial_dense_ev* solver, const char* solver_name, int argc, char** argv) {
  solver->ptr = new rokko::serial_dense_ev(std::string(solver_name));
  static_cast<rokko::serial_dense_ev*>(solver->ptr)->initialize(argc, argv);
}

void rokko_serial_dense_ev_construct_f(rokko_serial_dense_ev* solver, const char* solver_name) {
  int argc = 0;
  char** argv = nullptr;
  rokko_serial_dense_ev_construct(solver, solver_name, argc, argv);
}

void rokko_serial_dense_ev_destruct(rokko_serial_dense_ev* solver) {
  rokko::serial_dense_ev* ptr = static_cast<rokko::serial_dense_ev*>(solver->ptr);
  ptr->finalize();
  delete ptr;
  solver->ptr = nullptr;
}

// eigenvalues/eigenvectors, no parameters
struct rokko_parameters rokko_serial_dense_ev_diagonalize_eigen_matrix(struct rokko_serial_dense_ev solver,
									   struct rokko_eigen_matrix mat, struct rokko_eigen_vector eigvals,
									   struct rokko_eigen_matrix eigvecs) {
  struct rokko_parameters params_out;
  rokko_parameters_construct(&params_out);
  if (mat.major == rokko_matrix_col_major)
    *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::serial_dense_ev*>(solver.ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(mat.ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals.ptr),
								       *static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(eigvecs.ptr));
  else
    *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::serial_dense_ev*>(solver.ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(mat.ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals.ptr));
  return params_out;
}

// only eigenvalues, no parameters
struct rokko_parameters rokko_serial_dense_ev_diagonalize_eigvals(struct rokko_serial_dense_ev solver,
						   struct rokko_eigen_matrix mat, struct rokko_eigen_vector eigvals) {
  struct rokko_parameters params_out;
  rokko_parameters_construct(&params_out);
  if (mat.major == rokko_matrix_col_major)
    *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::serial_dense_ev*>(solver.ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(mat.ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals.ptr));
  else
    *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::serial_dense_ev*>(solver.ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(mat.ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals.ptr));
  return params_out;
}

// only eigenvalues, parameters
struct rokko_parameters rokko_serial_dense_ev_diagonalize_eigvals_params(struct rokko_serial_dense_ev solver,
    struct rokko_eigen_matrix mat, struct rokko_eigen_vector eigvals,
    struct rokko_parameters /* params */) {
  struct rokko_parameters params_out;
  rokko_parameters_construct(&params_out);
  if (mat.major == rokko_matrix_col_major)
    *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::serial_dense_ev*>(solver.ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(mat.ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals.ptr));
  else
    *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::serial_dense_ev*>(solver.ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(mat.ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals.ptr));
  return params_out;
}

// eigenvalues/eigenvectors, parameters
struct rokko_parameters rokko_serial_dense_ev_diagonalize_params(struct rokko_serial_dense_ev solver,
								 struct rokko_eigen_matrix mat, struct rokko_eigen_vector eigvals,
								 struct rokko_eigen_matrix eigvecs,
								 struct rokko_parameters params) {
  struct rokko_parameters params_out;
  rokko_parameters_construct(&params_out);
  if (mat.major == rokko_matrix_col_major)
    *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::serial_dense_ev*>(solver.ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(mat.ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals.ptr),
								       *static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(eigvecs.ptr),
								       *static_cast<rokko::parameters*>(params.ptr));
  else
    *static_cast<rokko::parameters*>(params_out.ptr) = static_cast<rokko::serial_dense_ev*>(solver.ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(mat.ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals.ptr),
								       *static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(eigvecs.ptr),
								       *static_cast<rokko::parameters*>(params.ptr));
  return params_out;
}

// For Fortran binding
// eigenvalues/eigenvectors, parameters
void rokko_serial_dense_ev_diagonalize_f(struct rokko_serial_dense_ev* solver,
					 struct rokko_eigen_matrix* mat, struct rokko_eigen_vector* eigvals,
					 struct rokko_eigen_matrix* eigvecs,
					 struct rokko_parameters* params, struct rokko_parameters* params_out) {
  rokko_parameters_construct(params_out);
  if (mat->major == rokko_matrix_col_major)
    *static_cast<rokko::parameters*>(params_out->ptr) = static_cast<rokko::serial_dense_ev*>(solver->ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(mat->ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								       *static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(eigvecs->ptr),
								       *static_cast<rokko::parameters*>(params->ptr));
  else
    *static_cast<rokko::parameters*>(params_out->ptr) = static_cast<rokko::serial_dense_ev*>(solver->ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(mat->ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								       *static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(eigvecs->ptr),
								       *static_cast<rokko::parameters*>(params->ptr));
}

void rokko_serial_dense_ev_diagonalize_no_params_out_f(struct rokko_serial_dense_ev* solver,
						       struct rokko_eigen_matrix* mat, struct rokko_eigen_vector* eigvals,
						       struct rokko_eigen_matrix* eigvecs,
						       struct rokko_parameters* params) {
  if (mat->major == rokko_matrix_col_major)
    static_cast<rokko::serial_dense_ev*>(solver->ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								   *static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(eigvecs->ptr),
								   *static_cast<rokko::parameters*>(params->ptr));
  else
    static_cast<rokko::serial_dense_ev*>(solver->ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								   *static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(eigvecs->ptr),
								   *static_cast<rokko::parameters*>(params->ptr));
}


void rokko_serial_dense_ev_diagonalize_no_params_inout_f(struct rokko_serial_dense_ev* solver,
							 struct rokko_eigen_matrix* mat, struct rokko_eigen_vector* eigvals,
							 struct rokko_eigen_matrix* eigvecs) {
  if (mat->major == rokko_matrix_col_major)
    static_cast<rokko::serial_dense_ev*>(solver->ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								   *static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(eigvecs->ptr));
  else
    static_cast<rokko::serial_dense_ev*>(solver->ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								   *static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(eigvecs->ptr));
}

void rokko_serial_dense_ev_diagonalize_eigvals_f(struct rokko_serial_dense_ev* solver,
						 struct rokko_eigen_matrix* mat, struct rokko_eigen_vector* eigvals,
						 struct rokko_parameters* params, struct rokko_parameters* params_out) {
  rokko_parameters_construct(params_out);
  if (mat->major == rokko_matrix_col_major)
    *static_cast<rokko::parameters*>(params_out->ptr) = static_cast<rokko::serial_dense_ev*>(solver->ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(mat->ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								       *static_cast<rokko::parameters*>(params->ptr));
  else
    *static_cast<rokko::parameters*>(params_out->ptr) = static_cast<rokko::serial_dense_ev*>(solver->ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(mat->ptr),
								       *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								       *static_cast<rokko::parameters*>(params->ptr));
}

void rokko_serial_dense_ev_diagonalize_eigvals_no_params_out_f(struct rokko_serial_dense_ev* solver,
							       struct rokko_eigen_matrix* mat, struct rokko_eigen_vector* eigvals,
							       struct rokko_parameters* params) {
  if (mat->major == rokko_matrix_col_major)
    static_cast<rokko::serial_dense_ev*>(solver->ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								   *static_cast<rokko::parameters*>(params->ptr));
  else
    static_cast<rokko::serial_dense_ev*>(solver->ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr),
								   *static_cast<rokko::parameters*>(params->ptr));
}

void rokko_serial_dense_ev_diagonalize_eigvals_no_params_inout_f(struct rokko_serial_dense_ev* solver,
								 struct rokko_eigen_matrix* mat, struct rokko_eigen_vector* eigvals) {
  if (mat->major == rokko_matrix_col_major)
    static_cast<rokko::serial_dense_ev*>(solver->ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr));
  else
    static_cast<rokko::serial_dense_ev*>(solver->ptr)->diagonalize(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(mat->ptr),
								   *static_cast<Eigen::VectorXd*>(eigvals->ptr));
}


int rokko_serial_dense_ev_num_solvers() {
  return rokko::serial_dense_ev::solvers().size();
}

char** rokko_serial_dense_ev_solvers() {
  std::vector<std::string> solvers = rokko::serial_dense_ev::solvers();
  char **solvers_c = (char**)malloc((size_t)(solvers.size() * sizeof(char*)));
  for (std::size_t i = 0; i < solvers.size(); ++i) {
    solvers_c[i] = copy_string(solvers[i]);
  }
  return solvers_c;
}

const char* rokko_serial_dense_ev_default_solver() {
  return rokko::serial_dense_ev::default_solver().c_str();
}
