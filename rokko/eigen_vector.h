/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_EIGEN_VECTOR_H
#define ROKKO_EIGEN_VECTOR_H

#ifdef __cplusplus
extern "C" {
#endif

struct rokko_eigen_vector {
  void* ptr;
};

void rokko_eigen_vector_construct(struct rokko_eigen_vector* vec, int dim1);
void rokko_eigen_vector_construct_array_size(struct rokko_eigen_vector* vector,
                                             int dim, double* ptr);
void rokko_eigen_vector_destruct(struct rokko_eigen_vector* vec);
int rokko_eigen_vector_get_dim(struct rokko_eigen_vector vec);
double* rokko_eigen_vector_get_array_pointer(struct rokko_eigen_vector vec);
void rokko_eigen_vector_set(struct rokko_eigen_vector vec, int i, double val);
void rokko_eigen_vector_set1(struct rokko_eigen_vector vec, int i, double val);
double rokko_eigen_vector_get(struct rokko_eigen_vector vec, int i);
double rokko_eigen_vector_get1(struct rokko_eigen_vector vec, int i);
void rokko_eigen_vector_print(struct rokko_eigen_vector vec);
void rokko_eigen_vector_generate_function(struct rokko_eigen_vector vec,
                                          double (*func)(int i));
void rokko_eigen_vector_generate_function_p(struct rokko_eigen_vector vec,
                                            double (*func)(const int* i));
void rokko_eigen_vector_generate_function1(struct rokko_eigen_vector vec,
                                            double (*func)(int i));
void rokko_eigen_vector_generate_function1_p(struct rokko_eigen_vector vec,
                                              double (*func)(const int* i));

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_EIGEN_VECTOR_H */
