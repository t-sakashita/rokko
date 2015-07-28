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

void rokko_serial_dense_solver_construct(rokko_serial_dense_solver* solver, const char* solver_name, int argc, char** argv) {
  solver->ptr = new rokko::serial_dense_solver(std::string(solver_name));
  static_cast<rokko::serial_dense_solver*>(solver->ptr)->initialize(argc, argv);
}

void rokko_serial_dense_solver_construct_f(rokko_serial_dense_solver* solver, const char* solver_name) {
  int argc = 0;
  char** argv;
  rokko_serial_dense_solver_construct(solver, solver_name, argc, argv);
}

void rokko_serial_dense_solver_destruct(rokko_serial_dense_solver* solver) {
  rokko::serial_dense_solver* ptr = static_cast<rokko::serial_dense_solver*>(solver->ptr);
  ptr->finalize();
  delete ptr;
}

// eigenvalues/eigenvectors, no parameters
void rokko_serial_dense_solver_diagonalize_localized_matrix(struct rokko_serial_dense_solver* solver,
							    struct rokko_localized_matrix* mat, struct rokko_localized_vector* eigvals,
							    struct rokko_localized_matrix* eigvecs) {
  if (mat->major == rokko_matrix_col_major)
    static_cast<rokko::serial_dense_solver*>(solver->ptr)->diagonalize(*static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(mat->ptr),
								       *static_cast<rokko::localized_vector<double>*>(eigvals->ptr),
								       *static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(eigvecs->ptr));
  else
    static_cast<rokko::serial_dense_solver*>(solver->ptr)->diagonalize(*static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(mat->ptr),
								       *static_cast<rokko::localized_vector<double>*>(eigvals->ptr));
}

// only eigenvalues, no parameters
void rokko_serial_dense_solver_diagonalize_eigvals(struct rokko_serial_dense_solver* solver,
						   struct rokko_localized_matrix* mat, struct rokko_localized_vector* eigvals) {
  if (mat->major == rokko_matrix_col_major)
    static_cast<rokko::serial_dense_solver*>(solver->ptr)->diagonalize(*static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(mat->ptr),
								       *static_cast<rokko::localized_vector<double>*>(eigvals->ptr));
  else
    static_cast<rokko::serial_dense_solver*>(solver->ptr)->diagonalize(*static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(mat->ptr),
								       *static_cast<rokko::localized_vector<double>*>(eigvals->ptr));
}

// only eigenvalues, parameters
void rokko_serial_dense_solver_diagonalize_eigvals_params(struct rokko_serial_dense_solver* solver,
							  struct rokko_localized_matrix* mat, struct rokko_localized_vector* eigvals,
							  struct rokko_parameters* params) {
  if (mat->major == rokko_matrix_col_major)
    static_cast<rokko::serial_dense_solver*>(solver->ptr)->diagonalize(*static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(mat->ptr),
								       *static_cast<rokko::localized_vector<double>*>(eigvals->ptr));
  else
    static_cast<rokko::serial_dense_solver*>(solver->ptr)->diagonalize(*static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(mat->ptr),
								       *static_cast<rokko::localized_vector<double>*>(eigvals->ptr));
}

// eigenvalues/eigenvectors, parameters
void rokko_serial_dense_solver_diagonalize_params(struct rokko_serial_dense_solver* solver,
						  struct rokko_localized_matrix* mat, struct rokko_localized_vector* eigvals,
						  struct rokko_localized_matrix* eigvecs,
						  struct rokko_parameters* params) {
  if (mat->major == rokko_matrix_col_major)
    static_cast<rokko::serial_dense_solver*>(solver->ptr)->diagonalize(*static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(mat->ptr),
								       *static_cast<rokko::localized_vector<double>*>(eigvals->ptr),
								       *static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(eigvecs->ptr),
								       *static_cast<rokko::parameters*>(params->ptr));
  else
    static_cast<rokko::serial_dense_solver*>(solver->ptr)->diagonalize(*static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(mat->ptr),
								       *static_cast<rokko::localized_vector<double>*>(eigvals->ptr),
								       *static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(eigvecs->ptr),
								       *static_cast<rokko::parameters*>(params->ptr));
  // std::list<std::string> keys = static_cast<rokko::parameters*>(params->ptr)->keys();
  //  for (std::list<std::string>::iterator i = keys.begin(); i != keys.end(); ++i) {
  //    std::cout << "keys=" << *i << std::endl;
  //  }
}

