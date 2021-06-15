/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2010-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

namespace rokko {

#include <algorithm>
#include <cctype>
#include <locale>
#include <string>

// trim from start (in place)
static inline void trim_left(std::string&& str) {
  str.erase(str.begin(), std::find_if_not(str.cbegin(), str.cend(), isspace));
}

// trim from end (in place)
inline void trim_right(std::string&& str) {
  str.erase(std::find_if_not(str.crbegin(), str.crend(), isspace).base(), str.end());
}

// trim from both ends (in place)
inline void trim(std::string&& str) {
  trim_left(std::move(str));
  trim_right(std::move(str));
}

// trim from start (copying)
inline std::string trim_left_copy(std::string&& str) {
  trim_left(std::move(str));
  return str;
}

// trim from end (copying)
inline std::string trim_right_copy(std::string&& str) {
  trim_right(std::move(str));
  return str;
}

// trim from both ends (copying)
inline std::string trim_copy(std::string&& str) {
  trim(std::move(str));
  return str;
}

} // namespace rokko
