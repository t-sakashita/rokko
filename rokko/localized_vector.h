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

#ifndef ROKKO_LOCALIZED_VECTOR_H
#define ROKKO_LOCALIZED_VECTOR_H

#ifdef __cplusplus
extern "C" {
#endif

struct rokko_localized_vector {
  void* ptr;
};

void rokko_localized_vector_construct(struct rokko_localized_vector* vec, int dim1);
void rokko_localized_vector_destruct(struct rokko_localized_vector* vec);
double rokko_localized_vector_get(struct rokko_localized_vector vec, int i);
double rokko_localized_vector_get_f(struct rokko_localized_vector vec, int i);
void rokko_localized_vector_print(struct rokko_localized_vector vec);

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_LOCALIZED_VECTOR_H */
