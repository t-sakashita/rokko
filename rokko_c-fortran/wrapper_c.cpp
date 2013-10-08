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
#include <rokko/collective.hpp>
#include "wrapper_c.h"

namespace rokko {

struct element_function_wrapper {
  element_function_wrapper(double (*f)(int*, int*)) { func_ptr = f; }
  double operator()(int i, int j) const {
    int is = i+1;
    int js = j+1;
    return func_ptr(&is, &js);
  }
  double (*func_ptr)(int*, int*);
};

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

void generate_distributed_matrix_function_row_major(void* mat, double (*func)(int i, int j)) {
  distributed_matrix<matrix_row_major>* mat_ = static_cast<distributed_matrix<matrix_row_major>*>(mat);
  mat_->generate(func);
}

void generate_distributed_matrix_function_col_major(void* mat, double (*func)(int i, int j)) {
  distributed_matrix<matrix_col_major>* mat_ = static_cast<distributed_matrix<matrix_col_major>*>(mat);
  mat_->generate(func);
}

void generate_distributed_matrix_function_row_major_fortran(void* mat, double (*func)(int* i, int* j)) {
  distributed_matrix<matrix_row_major>* mat_ = static_cast<distributed_matrix<matrix_row_major>*>(mat);
  element_function_wrapper wrapper(func);
  mat_->generate(wrapper);
}

void generate_distributed_matrix_function_col_major_fortran(void* mat, double (*func)(int* i, int* j)) {
  distributed_matrix<matrix_col_major>* mat_ = static_cast<distributed_matrix<matrix_col_major>*>(mat);
  element_function_wrapper wrapper(func);
  mat_->generate(wrapper);
}

// fortran array -> distributed_matrix
void generate_array_distributed_matrix_row_major_fortran(double* array, void* mat, int rows, int cols, int ld) {
  distributed_matrix<matrix_row_major>* mat_ = static_cast<distributed_matrix<matrix_row_major>*>(mat);
  for (int local_i=0; local_i<mat_->get_m_local(); ++local_i) {
    for (int local_j=0; local_j<mat_->get_n_local(); ++local_j) {
      int global_i = mat_->translate_l2g_row(local_i);
      int global_j = mat_->translate_l2g_col(local_j);
      mat_->set_local(local_i, local_j, array[global_i * ld + global_j]);
    }
  }
}

void generate_array_distributed_matrix_col_major_fortran(double* array, void* mat, int rows, int cols, int ld) {
  distributed_matrix<matrix_col_major>* mat_ = static_cast<distributed_matrix<matrix_col_major>*>(mat);
  for (int local_i=0; local_i<mat_->get_m_local(); ++local_i) {
    for (int local_j=0; local_j<mat_->get_n_local(); ++local_j) {
      int global_i = mat_->translate_l2g_row(local_i);
      int global_j = mat_->translate_l2g_col(local_j);
      mat_->set_local(local_i, local_j, array[global_i + global_j * ld]);
    }
  }
}

// distributed_matrix -> fortran array
void generate_distributed_matrix_array_row_major_fortran(void* mat, double* array, int rows, int cols, int ld) {
  distributed_matrix<matrix_row_major>* mat_ = static_cast<distributed_matrix<matrix_row_major>*>(mat);
  for (int local_i=0; local_i<mat_->get_m_local(); ++local_i) {
    for (int local_j=0; local_j<mat_->get_n_local(); ++local_j) {
      int global_i = mat_->translate_l2g_row(local_i);
      int global_j = mat_->translate_l2g_col(local_j);
      array[global_i * ld + global_j] = mat_->get_local(local_i, local_j);
    }
  }
}

void generate_distributed_matrix_array_col_major_fortran(void* mat, double* array, int rows, int cols, int ld) {
  distributed_matrix<matrix_col_major>* mat_ = static_cast<distributed_matrix<matrix_col_major>*>(mat);
  for (int local_i=0; local_i<mat_->get_m_local(); ++local_i) {
    for (int local_j=0; local_j<mat_->get_n_local(); ++local_j) {
      int global_i = mat_->translate_l2g_row(local_i);
      int global_j = mat_->translate_l2g_col(local_j);
      array[global_i + global_j * ld] = mat_->get_local(local_i, local_j);
    }
  }
}

void all_gather_fortran(void* mat, double* array) {
  distributed_matrix<rokko::matrix_col_major>* mat_ = static_cast<distributed_matrix<matrix_col_major>*>(mat);
  for(int root=0; root<mat_->get_nprocs(); ++root) {
    rokko::gather(*mat_, array, root);
  }
}

void set_distributed_matrix_local_col_major(void* mat, int i, int j, double val) {
  distributed_matrix<matrix_col_major>* mat_ = static_cast<distributed_matrix<matrix_col_major>*>(mat);
  mat_->set_local(i, j, val);
}

void set_distributed_matrix_local_row_major(void* mat, int i, int j, double val) {
  distributed_matrix<matrix_row_major>* mat_ = static_cast<distributed_matrix<matrix_row_major>*>(mat);
  mat_->set_local(i, j, val);
}

double get_distributed_matrix_local_col_major(void* mat, int i, int j) {
  distributed_matrix<matrix_col_major>* mat_ = static_cast<distributed_matrix<matrix_col_major>*>(mat);
  return mat_->get_local(i, j);
}

double get_distributed_matrix_local_row_major(void* mat, int i, int j) {
  distributed_matrix<matrix_row_major>* mat_ = static_cast<distributed_matrix<matrix_row_major>*>(mat);
  return mat_->get_local(i, j);
}

void print_distributed_matrix_col_major(void* mat) {
  distributed_matrix<matrix_col_major>* mat_ = static_cast<distributed_matrix<matrix_col_major>*>(mat);
  mat_->print();
}

void print_distributed_matrix_row_major(void* mat) {
  distributed_matrix<matrix_row_major>* mat_ = static_cast<distributed_matrix<matrix_row_major>*>(mat);
  mat_->print();
}

}
