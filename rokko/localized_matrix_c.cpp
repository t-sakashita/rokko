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

#include <rokko/localized_matrix.hpp>
#include <rokko/rokko_dense.h>

void rokko_localized_matrix_construct(rokko_localized_matrix* matrix, int dim1, int dim2,
				      int matrix_major) {
  if (matrix_major == rokko_matrix_col_major)
    matrix->ptr = new rokko::localized_matrix<double, rokko::matrix_col_major>(dim1, dim2);
  else
    matrix->ptr = new rokko::localized_matrix<double, rokko::matrix_row_major>(dim1, dim2);
  matrix->major = matrix_major;
}

void rokko_localized_matrix_destruct(rokko_localized_matrix* matrix) {
  if (matrix->major == rokko_matrix_col_major)
    delete static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(matrix->ptr);
  else
    delete static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(matrix->ptr);
}

void rokko_localized_matrix_print(rokko_localized_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->print();
  else
    static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->print();
}


void rokko_localized_matrix_generate_function(struct rokko_localized_matrix* matrix,
					      double (*func)(int i, int j)) {
  if (matrix->major == rokko_matrix_col_major)
    static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(matrix->ptr)->generate(func);
  else
    static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(matrix->ptr)->generate(func);
}

void rokko_localized_matrix_set_local(rokko_localized_matrix* matrix, int local_i, int local_j, double value) {
  if (matrix->major == rokko_matrix_col_major)
    static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(matrix->ptr)->set_local(local_i,local_j,value);
  else
    static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(matrix->ptr)->set_local(local_i,local_j,value);
}

double rokko_localized_matrix_get_local(rokko_localized_matrix matrix, int local_i, int local_j) { 
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_local(local_i,local_j);
  else
    return static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_local(local_i,local_j);
}

void rokko_localized_matrix_set_global(rokko_localized_matrix* matrix, int global_i, int global_j, double value) {
  if (matrix->major == rokko_matrix_col_major)
    static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(matrix->ptr)->set_global(global_i,global_j,value);
  else
    static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(matrix->ptr)->set_global(global_i,global_j,value);
}

double rokko_localized_matrix_get_global(rokko_localized_matrix matrix, int global_i, int global_j) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_global(global_i,global_j);
  else
    return static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_global(global_i,global_j);
}

int rokko_localized_matrix_get_m_local(struct rokko_localized_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_m_local();
  else
    return static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_m_local();
}

int rokko_localized_matrix_get_n_local(struct rokko_localized_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_n_local();
  else
    return static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_n_local();
}

int rokko_localized_matrix_get_m_global(struct rokko_localized_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_m_global();
  else
    return static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_m_global();
}

int rokko_localized_matrix_get_n_global(struct rokko_localized_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::localized_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_n_global();
  else
    return static_cast<rokko::localized_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_n_global();
}

