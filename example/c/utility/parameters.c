/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/parameters.h>
#include <stdio.h>

int main(/* int argc , char** argv */) {
  struct rokko_parameters params;
  rokko_parameters_construct(&params);
  char p = 'V';
  // set some parameters
  rokko_parameters_set_double(params, "T", 0.95);
  rokko_parameters_set_double(params, "Pi", 3.14);
  rokko_parameters_set_string(params, "A", "test");
  rokko_parameters_set_char(params, "C", 'A');
  rokko_parameters_set_char(params, "C", p);

  // get double
  rokko_parameters_get_double(params, "T");

  // get value as string
  printf("value=%s\n", rokko_parameters_get_string(params, "A"));

  // get char
  rokko_parameters_get_char(params, "C");

  // is "T" defined?
  rokko_parameters_defined(params, "T");

  // is "S" defined?
  rokko_parameters_defined(params, "S");

  // clear paramter "Pi"
  rokko_parameters_clear(params, "Pi");

  // output list of parameters as string
  int num_keys = rokko_parameters_size(params);
  char** keys = rokko_parameters_keys(params);
  int i;
  for (i = 0; i < num_keys; ++i) {
    printf("%s = %s\n", keys[i], rokko_parameters_get_string(params, keys[i]));
  }

  rokko_parameters_destruct(&params);
}
