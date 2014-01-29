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
#include <rokko/distributed_matrix.hpp>
#include <rokko/collective.hpp>
#include <rokko/rokko.h>

void rokko_distributed_matrix_construct(rokko_distributed_matrix* matrix, int dim1, int dim2,
                                        rokko_grid grid, rokko_solver solver, int matrix_major) {
  if (matrix_major == rokko_matrix_col_major)
    matrix->ptr = new rokko::distributed_matrix<rokko::matrix_col_major>(dim1, dim2,
      *static_cast<rokko::grid*>(grid.ptr), *static_cast<rokko::solver*>(solver.ptr));
  else
    matrix->ptr = new rokko::distributed_matrix<rokko::matrix_row_major>(dim1, dim2,
      *static_cast<rokko::grid*>(grid.ptr), *static_cast<rokko::solver*>(solver.ptr));
  matrix->major = matrix_major;
}
void rokko_distributed_matrix_destruct(rokko_distributed_matrix* matrix) {
  if (matrix->major == rokko_matrix_col_major)
    delete static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix->ptr);
  else
    delete static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix->ptr);
}

void rokko_distributed_matrix_generate_function(struct rokko_distributed_matrix* matrix,
  double (*func)(int i, int j)) {
  if (matrix->major == rokko_matrix_col_major)
    static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix->ptr)->generate(func);
  else
    static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix->ptr)->generate(func);
}

void rokko_distributed_matrix_print(rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix.ptr)->print();
  else
    static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix.ptr)->print();
}

void rokko_distributed_matrix_set_local(rokko_distributed_matrix* matrix, int local_i, int local_j, double value) {
    if (matrix->major == rokko_matrix_col_major)
      static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix->ptr)->set_local(local_i,local_j,value);
    else
      static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix->ptr)->set_local(local_i,local_j,value);
  }

double rokko_distributed_matrix_get_local(rokko_distributed_matrix matrix, int local_i, int local_j)  { 
    if (matrix.major == rokko_matrix_col_major)
      return static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix.ptr)->get_local(local_i,local_j);
    else
      return static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix.ptr)->get_local(local_i,local_j);
  }

void rokko_distributed_matrix_set_global(rokko_distributed_matrix* matrix, int global_i, int global_j, double value) {
    if (matrix->major == rokko_matrix_col_major)
      static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix->ptr)->set_global(global_i,global_j,value);
    else
      static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix->ptr)->set_global(global_i,global_j,value);
  }

double rokko_distributed_matrix_get_global(rokko_distributed_matrix matrix, int global_i, int global_j) {
    if (matrix.major == rokko_matrix_col_major)
      return static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix.ptr)->get_global(global_i,global_j);
    else
      return static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix.ptr)->get_global(global_i,global_j);
  }

int rokko_distributed_matrix_get_m_local(struct rokko_distributed_matrix matrix){
    if (matrix.major == rokko_matrix_col_major)
      return static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix.ptr)->get_m_local();
    else
      return static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix.ptr)->get_m_local();
  }
int rokko_distributed_matrix_get_n_local(struct rokko_distributed_matrix matrix){
    if (matrix.major == rokko_matrix_col_major)
      return static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix.ptr)->get_n_local();
    else
      return static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix.ptr)->get_n_local();
  }
int rokko_distributed_matrix_get_m_global(struct rokko_distributed_matrix matrix){
    if (matrix.major == rokko_matrix_col_major)
      return static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix.ptr)->get_m_global();
    else
      return static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix.ptr)->get_m_global();
  }
int rokko_distributed_matrix_get_n_global(struct rokko_distributed_matrix matrix){
    if (matrix.major == rokko_matrix_col_major)
      return static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix.ptr)->get_n_global();
    else
      return static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix.ptr)->get_n_global();
  }

int rokko_distributed_matrix_get_nprocs(struct rokko_distributed_matrix matrix){
    if (matrix.major == rokko_matrix_col_major)
      return static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix.ptr)->get_nprocs();
    else
      return static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix.ptr)->get_nprocs();
  }

int rokko_distributed_matrix_get_myrank(struct rokko_distributed_matrix matrix){
    if (matrix.major == rokko_matrix_col_major)
      return static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix.ptr)->get_myrank();
    else
      return static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix.ptr)->get_myrank();
  }

int rokko_distributed_matrix_translate_l2g_row(struct rokko_distributed_matrix matrix, int local_i){
    if (matrix.major == rokko_matrix_col_major)
      return static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix.ptr)->translate_l2g_row(local_i);
    else
      return static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix.ptr)->translate_l2g_row(local_i);
  }
int rokko_distributed_matrix_translate_l2g_col(struct rokko_distributed_matrix matrix, int local_j){
    if (matrix.major == rokko_matrix_col_major)
      return static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix.ptr)->translate_l2g_col(local_j);
    else
      return static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix.ptr)->translate_l2g_col(local_j);
  }
int rokko_distributed_matrix_translate_g2l_row(struct rokko_distributed_matrix matrix, int global_i){
    if (matrix.major == rokko_matrix_col_major)
      return static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix.ptr)->translate_g2l_row(global_i);
    else
      return static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix.ptr)->translate_g2l_row(global_i);
  }
int rokko_distributed_matrix_translate_g2l_col(struct rokko_distributed_matrix matrix, int global_j){
    if (matrix.major == rokko_matrix_col_major)
      return static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix.ptr)->translate_g2l_col(global_j);
    else
      return static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix.ptr)->translate_g2l_col(global_j);
  }

int rokko_gather(struct rokko_distributed_matrix* matrix, double* array, int root){
  if (matrix->major == rokko_matrix_col_major){
    rokko::distributed_matrix<rokko::matrix_col_major>* ptr_ = static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix->ptr);
    return rokko::gather(*ptr_, array, root);
  }
  else{
    rokko::distributed_matrix<rokko::matrix_row_major>* ptr_ = static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix->ptr);
    return rokko::gather(*ptr_, array, root);
  }
}

int rokko_scatter(double* global_array, struct rokko_distributed_matrix* matrix, int root){
  if (matrix->major == rokko_matrix_col_major){
    rokko::distributed_matrix<rokko::matrix_col_major>* ptr_ = static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix->ptr);
    return rokko::scatter(global_array, *ptr_, root);
  }
  else{
    rokko::distributed_matrix<rokko::matrix_row_major>* ptr_ = static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix->ptr);
    return rokko::scatter(global_array, *ptr_, root);
  }
}

void rokko_all_gather(struct rokko_distributed_matrix* matrix, double* array){
  if (matrix->major == rokko_matrix_col_major){
    rokko::distributed_matrix<rokko::matrix_col_major>* ptr_ = static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix->ptr);
    for(int root = 0; root < ptr_->get_nprocs(); ++root) {
      rokko::gather(*ptr_, array, root);
    }
  }
  else{
    rokko::distributed_matrix<rokko::matrix_row_major>* ptr_ = static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix->ptr);
    for(int root = 0; root < ptr_->get_nprocs(); ++root) {
      rokko::gather(*ptr_, array, root);
    }
  }
}
