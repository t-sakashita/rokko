/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/distributed_matrix.h>
#include <rokko/distributed_matrix.hpp>

void rokko_distributed_matrix_construct(struct rokko_distributed_matrix* matrix, struct rokko_mapping_bc map) {
  if (map.major == rokko_matrix_col_major)
    matrix->ptr = new rokko::distributed_matrix<double, rokko::matrix_col_major>(*static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr));
  else
    matrix->ptr = new rokko::distributed_matrix<double, rokko::matrix_row_major>(*static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr));
  matrix->major = map.major;
}

void rokko_distributed_matrix_construct_array(struct rokko_distributed_matrix* matrix, struct rokko_mapping_bc map, double *array) {
  if (map.major == rokko_matrix_col_major)
    matrix->ptr = new rokko::distributed_matrix<double, rokko::matrix_col_major>(*static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr), array);
  else
    matrix->ptr = new rokko::distributed_matrix<double, rokko::matrix_row_major>(*static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr), array);
  matrix->major = map.major;
}

void rokko_distributed_matrix_construct_array_sizes(struct rokko_distributed_matrix* matrix, struct rokko_mapping_bc map, int dim1, int dim2, double *array) {
  if (map.major == rokko_matrix_col_major)
    matrix->ptr = new rokko::distributed_matrix<double, rokko::matrix_col_major>(*static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr), array);
  else
    matrix->ptr = new rokko::distributed_matrix<double, rokko::matrix_row_major>(*static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr), array);
  matrix->major = map.major;
}

void rokko_distributed_matrix_destruct(struct rokko_distributed_matrix* matrix) {
  if (matrix->major == rokko_matrix_col_major)
    delete static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix->ptr);
  else
    delete static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix->ptr);
  matrix->ptr = nullptr;
}

void rokko_distributed_matrix_generate_function(struct rokko_distributed_matrix matrix,
						double (*func)(int i, int j)) {
  if (matrix.major == rokko_matrix_col_major)
    static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->generate(func);
  else
    static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->generate(func);
}

void rokko_distributed_matrix_generate_function_p(struct rokko_distributed_matrix matrix,
						double (*func)(const int* i, const int* j)) {
  const auto g = [&func](int i, int j) { return func(&i, &j); };

  if (matrix.major == rokko_matrix_col_major)
    static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->generate(g);
  else
    static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->generate(g);
}

void rokko_distributed_matrix_generate_function1(struct rokko_distributed_matrix matrix,
						double (*func)(int i, int j)) {
  auto const g = [&func](int i, int j) { return func(i+1, j+1); };

  if (matrix.major == rokko_matrix_col_major)
    static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->generate(g);
  else
    static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->generate(g);
}

void rokko_distributed_matrix_generate_function1_p(struct rokko_distributed_matrix matrix,
						double (*func)(const int* i, const int* j)) {
  const auto g = [&func](int i, int j) {
    const int i1 = i+1, j1 = j+1;
    return func(&i1, &j1);
  };
  
  if (matrix.major == rokko_matrix_col_major)
    static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->generate(g);
  else
    static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->generate(g);
}

void rokko_distributed_matrix_print(rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->print();
  else
    static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->print();
}

void rokko_distributed_matrix_set_local(rokko_distributed_matrix matrix, int local_i, int local_j, double value) {
  if (matrix.major == rokko_matrix_col_major)
    static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->set_local(local_i,local_j,value);
  else
    static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->set_local(local_i,local_j,value);
}

double rokko_distributed_matrix_get_local(rokko_distributed_matrix matrix, int local_i, int local_j) { 
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_local(local_i,local_j);
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_local(local_i,local_j);
}

void rokko_distributed_matrix_set_global(rokko_distributed_matrix matrix, int global_i, int global_j, double value) {
  if (matrix.major == rokko_matrix_col_major)
    static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->set_global(global_i,global_j,value);
  else
    static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->set_global(global_i,global_j,value);
}

double rokko_distributed_matrix_get_global(rokko_distributed_matrix matrix, int global_i, int global_j) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_global(global_i,global_j);
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_global(global_i,global_j);
}

void rokko_distributed_matrix_set_local1(rokko_distributed_matrix matrix, int local_i, int local_j, double value) {
  if (matrix.major == rokko_matrix_col_major)
    static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->set_local(local_i-1,local_j-1,value);
  else
    static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->set_local(local_i-1,local_j-1,value);
}

double rokko_distributed_matrix_get_local1(rokko_distributed_matrix matrix, int local_i, int local_j) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_local(local_i-1,local_j-1);
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_local(local_i-1,local_j-1);
}

void rokko_distributed_matrix_set_global1(rokko_distributed_matrix matrix, int global_i, int global_j, double value) {
  if (matrix.major == rokko_matrix_col_major)
    static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->set_global(global_i-1,global_j-1,value);
  else
    static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->set_global(global_i-1,global_j-1,value);
}

double rokko_distributed_matrix_get_global1(rokko_distributed_matrix matrix, int global_i, int global_j) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_global(global_i-1,global_j-1);
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_global(global_i-1,global_j-1);
}

int rokko_distributed_matrix_get_mb(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_mb();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_mb();
}

int rokko_distributed_matrix_get_nb(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_nb();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_nb();
}

int rokko_distributed_matrix_get_m_local(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_m_local();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_m_local();
}

int rokko_distributed_matrix_get_n_local(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_n_local();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_n_local();
}

int rokko_distributed_matrix_get_m_global(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_m_global();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_m_global();
}

int rokko_distributed_matrix_get_n_global(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_n_global();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_n_global();
}

int rokko_distributed_matrix_get_lld(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_lld();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_lld();
}

int rokko_distributed_matrix_get_m_size(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_m_size();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_m_size();
}

int rokko_distributed_matrix_get_n_size(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_n_size();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_n_size();
}

int rokko_distributed_matrix_get_nprocs(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_nprocs();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_nprocs();
}

int rokko_distributed_matrix_get_myrank(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_myrank();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_myrank();
}

int rokko_distributed_matrix_get_nprow(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_nprow();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_nprow();
}

int rokko_distributed_matrix_get_npcol(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_npcol();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_npcol();
}

int rokko_distributed_matrix_get_myrow(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_myrow();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_myrow();
}

int rokko_distributed_matrix_get_mycol(struct rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_mycol();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_mycol();
}

int rokko_distributed_matrix_translate_l2g_row(struct rokko_distributed_matrix matrix, int local_i) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->translate_l2g_row(local_i);
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->translate_l2g_row(local_i);
}

int rokko_distributed_matrix_translate_l2g_col(struct rokko_distributed_matrix matrix, int local_j) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->translate_l2g_col(local_j);
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->translate_l2g_col(local_j);
}

int rokko_distributed_matrix_translate_g2l_row(struct rokko_distributed_matrix matrix, int global_i) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->translate_g2l_row(global_i);
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->translate_g2l_row(global_i);
}

int rokko_distributed_matrix_translate_g2l_col(struct rokko_distributed_matrix matrix, int global_j) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->translate_g2l_col(global_j);
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->translate_g2l_col(global_j);
}

// offset by one for Fortran
int rokko_distributed_matrix_translate_l2g_row1(struct rokko_distributed_matrix matrix, int local_i) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->translate_l2g_row(local_i - 1) + 1;
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->translate_l2g_row(local_i - 1) + 1;
}

int rokko_distributed_matrix_translate_l2g_col1(struct rokko_distributed_matrix matrix, int local_j) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->translate_l2g_col(local_j - 1) + 1;
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->translate_l2g_col(local_j - 1) + 1;
}

int rokko_distributed_matrix_translate_g2l_row1(struct rokko_distributed_matrix matrix, int global_i) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->translate_g2l_row(global_i - 1) + 1;
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->translate_g2l_row(global_i - 1) + 1;
}

int rokko_distributed_matrix_translate_g2l_col1(struct rokko_distributed_matrix matrix, int global_j) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->translate_g2l_col(global_j - 1) + 1;
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->translate_g2l_col(global_j - 1) + 1;
}

bool rokko_distributed_matrix_is_row_major(struct rokko_distributed_matrix matrix) {
  return matrix.major == rokko_matrix_row_major;
}

bool rokko_distributed_matrix_is_col_major(struct rokko_distributed_matrix matrix) {
  return matrix.major == rokko_matrix_col_major;
}

double* rokko_distributed_matrix_get_array_pointer(struct rokko_distributed_matrix matrix) { 
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr)->get_array_pointer();
  else
    return static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr)->get_array_pointer();
}


