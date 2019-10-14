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

#include <rokko/solver.hpp>
#include <rokko/mapping_bc.hpp>
#include <rokko/dense.h>
#include <rokko/mapping_bc.h>

void rokko_mapping_bc_construct_block_size(struct rokko_mapping_bc* map, int global_dim, int block_size, struct rokko_grid grid) {
  // Fix me: All parallel dense solvers use matrix_col_major. So we omit its checking, so far...
  map->ptr = new rokko::mapping_bc<rokko::matrix_col_major>(global_dim, block_size, *static_cast<rokko::grid*>(grid.ptr));
  map->major = rokko_matrix_col_major;
}

void rokko_mapping_bc_destruct(struct rokko_mapping_bc* map) {
  delete static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map->ptr);
  map->ptr = nullptr;
}

int rokko_mapping_bc_get_m_local(struct rokko_mapping_bc matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(matrix.ptr)->get_m_local();
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(matrix.ptr)->get_m_local();
}

int rokko_mapping_bc_get_n_local(struct rokko_mapping_bc matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(matrix.ptr)->get_n_local();
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(matrix.ptr)->get_n_local();
}

int rokko_mapping_bc_get_m_global(struct rokko_mapping_bc matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(matrix.ptr)->get_m_global();
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(matrix.ptr)->get_m_global();
}

int rokko_mapping_bc_get_n_global(struct rokko_mapping_bc matrix) {
  if (matrix.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(matrix.ptr)->get_n_global();
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(matrix.ptr)->get_n_global();
}

int rokko_mapping_bc_get_m_size(struct rokko_mapping_bc map) {
  if (map.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr)->get_m_size();
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr)->get_m_size();
}

int rokko_mapping_bc_get_n_size(struct rokko_mapping_bc map) {
  if (map.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr)->get_n_size();
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr)->get_n_size();
}

int rokko_mapping_bc_get_length_array(struct rokko_mapping_bc map) {
  if (map.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr)->get_length_array();
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr)->get_length_array();
}

int rokko_mapping_bc_get_nprocs(struct rokko_mapping_bc map) {
  if (map.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr)->get_nprocs();
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr)->get_nprocs();
}

int rokko_mapping_bc_get_myrank(struct rokko_mapping_bc map) {
  if (map.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr)->get_myrank();
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr)->get_myrank();
}

int rokko_mapping_bc_get_nprow(struct rokko_mapping_bc map) {
  if (map.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr)->get_nprow();
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr)->get_nprow();
}

int rokko_mapping_bc_get_npcol(struct rokko_mapping_bc map) {
  if (map.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr)->get_npcol();
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr)->get_npcol();
}

int rokko_mapping_bc_get_myrow(struct rokko_mapping_bc map) {
  if (map.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr)->get_myrow();
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr)->get_myrow();
}

int rokko_mapping_bc_get_mycol(struct rokko_mapping_bc map) {
  if (map.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr)->get_mycol();
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr)->get_mycol();
}

int rokko_mapping_bc_translate_l2g_row(struct rokko_mapping_bc map, int local_i) {
  if (map.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr)->translate_l2g_row(local_i);
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr)->translate_l2g_row(local_i);
}

int rokko_mapping_bc_translate_l2g_col(struct rokko_mapping_bc map, int local_j) {
  if (map.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr)->translate_l2g_col(local_j);
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr)->translate_l2g_col(local_j);
}

int rokko_mapping_bc_translate_g2l_row(struct rokko_mapping_bc map, int global_i) {
  if (map.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr)->translate_g2l_row(global_i);
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr)->translate_g2l_row(global_i);
}

int rokko_mapping_bc_translate_g2l_col(struct rokko_mapping_bc map, int global_j) {
  if (map.major == rokko_matrix_col_major)
    return static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map.ptr)->translate_g2l_col(global_j);
  else
    return static_cast<rokko::mapping_bc<rokko::matrix_row_major>*>(map.ptr)->translate_g2l_col(global_j);
}

bool rokko_mapping_bc_is_row_major(struct rokko_mapping_bc map) {
  if (map.major == rokko_matrix_row_major)
    return true;
  else
    return false;
}

bool rokko_mapping_bc_is_col_major(struct rokko_mapping_bc map) {
  if (map.major == rokko_matrix_col_major)
    return true;
  else
    return false;
}
