/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_EIGEN_MATRIX_H
#define ROKKO_EIGEN_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

struct rokko_eigen_matrix {
  void* ptr;
  int major;
};

void rokko_eigen_matrix_construct(struct rokko_eigen_matrix* matrix, int dim1, int dim2,
  int matrix_major);
void rokko_eigen_matrix_destruct(struct rokko_eigen_matrix* matrix);
void rokko_eigen_matrix_print(struct rokko_eigen_matrix matrix);
void rokko_eigen_matrix_generate_function(struct rokko_eigen_matrix matrix,
  double (*func)(int i, int j));
void rokko_eigen_matrix_set_local(struct rokko_eigen_matrix matrix,
  int local_i, int local_j, double value);
double rokko_eigen_matrix_get_local(struct rokko_eigen_matrix matrix,
  int local_i, int local_j);
void rokko_eigen_matrix_set_global(struct rokko_eigen_matrix matrix,
  int global_i, int global_j, double value);
double rokko_eigen_matrix_get_global(struct rokko_eigen_matrix matrix,
  int global_i, int global_j);
int rokko_eigen_matrix_get_m_local(struct rokko_eigen_matrix matrix);
int rokko_eigen_matrix_get_n_local(struct rokko_eigen_matrix matrix);
int rokko_eigen_matrix_get_m_global(struct rokko_eigen_matrix matrix);
int rokko_eigen_matrix_get_n_global(struct rokko_eigen_matrix matrix);
double* rokko_eigen_matrix_get_pointer(struct rokko_eigen_matrix matrix);

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_EIGEN_MATRIX_H */
