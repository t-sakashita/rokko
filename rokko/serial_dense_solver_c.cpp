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

#include <rokko/solver.hpp>
#include <rokko/rokko_dense.h>

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

void rokko_serial_dense_solver_diagonalize_localized_matrix(struct rokko_serial_dense_solver* solver,
							    struct rokko_localized_matrix* mat, struct rokko_localized_vector* eigvals,
							    struct rokko_localized_matrix* eigvecs) {
  if (mat->major == rokko_matrix_col_major)
    static_cast<rokko::serial_dense_solver*>(solver->ptr)->diagonalize(*static_cast<rokko::localized_matrix<rokko::matrix_col_major>*>(mat->ptr),
								       *static_cast<rokko::localized_vector*>(eigvals->ptr),
								       *static_cast<rokko::localized_matrix<rokko::matrix_col_major>*>(eigvecs->ptr));
  else
    static_cast<rokko::serial_dense_solver*>(solver->ptr)->diagonalize(*static_cast<rokko::localized_matrix<rokko::matrix_row_major>*>(mat->ptr),
								       *static_cast<rokko::localized_vector*>(eigvals->ptr),
								       *static_cast<rokko::localized_matrix<rokko::matrix_row_major>*>(eigvecs->ptr));
}
