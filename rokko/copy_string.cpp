/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <string>

char* copy_string(std::string const& str) {
  char* p = (char*) malloc( sizeof(char) * (str.size() + 1) );
  str.copy(p, str.size(), 0);
  p[str.size()] = '\0';
  return p;
}
