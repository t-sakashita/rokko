/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>,
*                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>
#include "wrapper_c.h"

namespace rokko {

void* initialize_distributed_matrix_row_major(int dim1, int dim2, void* g, void* solver_ ) {
  grid* g_ = static_cast<grid*>(g);
  solver* solver__ =static_cast<solver*>(solver_);
  return static_cast<void*>(new distributed_matrix<matrix_row_major>(dim1, dim2, *g_, *solver__));
}

void* initialize_distributed_matrix_col_major(int dim1, int dim2, void* g, void* solver_ ) {
  grid* g_ = static_cast<grid*>(g);
  solver* solver__ =static_cast<solver*>(solver_);
  return static_cast<void*>(new distributed_matrix<matrix_col_major>(dim1, dim2, *g_, *solver__));
}

void delete_distributed_matrix_col_major(void*  mat) {
  distributed_matrix<matrix_col_major>* mat_ = static_cast<distributed_matrix<matrix_col_major>* >(mat);
  delete mat_;
}

void delete_distributed_matrix_row_major(void*  mat) {
  distributed_matrix<matrix_row_major>* mat_ = static_cast<distributed_matrix<matrix_row_major>* >(mat);
  delete mat_;
}

void* initialize_localized_vector(int dim) {
  return static_cast<void*>(new localized_vector(dim));
}

void delete_localized_vector(void* w ) {
  localized_vector* w_ = static_cast<localized_vector*>(w);
  delete w_;
}

double localized_vector_get_element(void* w, int i) {
  localized_vector* w_ = static_cast<localized_vector*>(w);
  return (*w_)(i);
}

void* initialize_grid_col_major(MPI_Comm comm) {
  return static_cast<void*> (new grid(comm, grid_col_major));
}

void* initialize_grid_row_major(MPI_Comm comm) {
  grid_row_major_t grid_row_major;
  return static_cast<void*> (new grid(comm, grid_row_major));
}

int grid_get_myrank(void* g) {
  grid* g_= static_cast<grid*>(g);
  return g_->get_myrank();
}

int grid_get_nprocs(void* g) {
  grid* g_= static_cast<grid*>(g);
  return g_->get_nprocs();
}

void delete_grid(void* g) {
  grid* g_= static_cast<grid*>(g);
  delete g_;
}

void* initialize_solver(char* solver_name, int argc, char *argv[]) {
  std::string solver_name_in(solver_name);
  solver* solver_ = new solver(solver_name_in);
  solver_->initialize(argc,argv);
  return static_cast<void*>(solver_);
}

void delete_solver(void* solver_) {
  solver* solver__ = static_cast<solver*>(solver_);
  solver__ -> finalize();
  delete solver__;
}

void solver_diagonalize_matrix_col_major(void* solver_ ,void* mat, void* w, void* Z, void* timer) {
  solver* solver__ = static_cast<solver*>(solver_);
  distributed_matrix<matrix_col_major>* mat_ = static_cast<distributed_matrix<matrix_col_major>*> (mat);
  localized_vector* w_ = static_cast<localized_vector*>(w);
  distributed_matrix<matrix_col_major>* Z_ = static_cast<distributed_matrix<matrix_col_major>*> (Z);
  timer_dumb* timer_ = static_cast<timer_dumb*>(timer);
  solver__ -> diagonalize(*mat_, *w_, *Z_, *timer_);
}

void solver_diagonalize_matrix_row_major(void* solver_ ,void* mat, void* w, void* Z, void* timer) {
  solver* solver__ = static_cast<solver*>(solver_);
  distributed_matrix<matrix_row_major>* mat_ = static_cast<distributed_matrix<matrix_row_major>*> (mat);
  localized_vector* w_ = static_cast<localized_vector*>(w);
  distributed_matrix<matrix_row_major>* Z_ = static_cast<distributed_matrix<matrix_row_major>*> (Z);
  timer_dumb* timer_ = static_cast<timer_dumb*>(timer);

  solver__ -> diagonalize(*mat_, *w_, *Z_, *timer_);
}

void generate_distributed_matrix_col_major(void* mat, double (*func)(int i, int j)) {
  distributed_matrix<matrix_col_major>* mat_ = static_cast<distributed_matrix<matrix_col_major>*>(mat);
  std::cout << "m_gglobal=" << mat_->get_m_global() << "n_gglobal=" << mat_->get_n_global() << std::endl;
  mat_->generate(func);
}

}
