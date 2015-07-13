/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_PARAMETERS_H
#define ROKKO_PARAMETERS_H

#include <rokko/config.h>

#ifdef __cplusplus
extern "C" {
#endif

struct rokko_parameters {
  void* ptr;
};

void rokko_parameters_construct(struct rokko_parameters* params);
void rokko_parameters_destruct(struct rokko_parameters* params);

void rokko_parameters_clear_all(struct rokko_parameters* params);
void rokko_parameters_clear(struct rokko_parameters* params, const char* key);

void rokko_parameters_set_int(struct rokko_parameters* params, const char* key, int value);

void rokko_parameters_set_double(struct rokko_parameters* params, const char* key, double value);

void rokko_parameters_set_char(struct rokko_parameters* params, const char* key, char value);

void rokko_parameters_set_string(struct rokko_parameters* params, const char* key, const char* value);

int rokko_parameters_get_int(struct rokko_parameters* params, const char* key);

double rokko_parameters_get_double(struct rokko_parameters* params, const char* key);

char rokko_parameters_get_char(struct rokko_parameters* params, const char* key);

char* rokko_parameters_get_string(struct rokko_parameters* params, const char* key);

int rokko_parameters_defined(struct rokko_parameters* params, const char* key);

int rokko_parameters_size(struct rokko_parameters* params);

char** rokko_parameters_keys(struct rokko_parameters* params);

#ifdef __cplusplus
}
#endif

#endif /* ROKKO_PARAMETERS_H */
